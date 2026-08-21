# test_persistence_flag.py

import json

import asdf
import pytest
import numpy as np

from astropy import coordinates, units
from astropy.io import fits
from astropy.modeling import models
from gwcs import coordinate_frames, wcs
from psfsim.polychrom import PolychromaticPSF
from roman_datamodels.dqflags import pixel
from scipy.ndimage import shift

from roman_hlis_l2_driver.persistence.persistence_flag import previous_obsid, get_prev_obs
from roman_hlis_l2_driver.persistence.persistence import run


def make_wcs(ra=9.55, dec=-44.1, lonpole=180.0):
    """ Creates a GWCS for mock asdf files.
    """
    shiftscale = (models.Shift(-2043.5) | models.Scale(0.11/3600.0))
    det2sky = ( ((shiftscale | models.Scale(-1)) & shiftscale)
               | models.Pix2Sky_TAN() | models.RotateNative2Celestial(ra,dec,lonpole) )
    
    return wcs.WCS(
        [
            (
                coordinate_frames.Frame2D(
                    name="detector", axes_names=("x", "y"), unit=(units.pix, units.pix)
                ),
                det2sky,
            ),
            (
                coordinate_frames.CelestialFrame(
                    reference_frame=coordinates.ICRS(), name="icrs", unit=(units.deg, units.deg)
                ),
                None,
            ),
        ]
    )

def inject_star(image,wcsobj,psf,ra,dec,flux,ov=1,nbdy=256,n=121):
    """ Inject a star into a mock image
    """
    # where is this?
    x, y = wcsobj.invert(ra_s, dec_s)
    xc = int(np.round(x))
    yc = int(np.round(y))

    psf_shifted = shift(psf, (ov * (x - xc), ov * (y - yc)), order=1, mode="constant")
    stamp = np.sum(psf_shifted.reshape((n, ov, n, ov)), axis=(1, 3) )
    #update the image object that was passed
    image[
        yc + nbdy - n // 2 : yc + nbdy + n // 2 + 1, xc + nbdy - n // 2 : xc + nbdy + n // 2 + 1
    ] += flux_s * stamp
    
    return yc, xc


@pytest.fixture
def obsfile(tmp_path):
    """
    """
    obsfile = tmp_path/"Roman_WAS_obseq_test.fits"
    n_obs = 10
    exptime = 100.0
    cadence_sec = 180.0
    cadence_day = cadence_sec / 86400.0
    
    dates = (62000.0 + np.arrange(n_obs,dtype=np.float64) * cadence_day)
    columns = [
        fits.Column(
            name="date",
            format="D",
            array=dates,
        ),
        fits.Column(
            name="exptime",
            format="D",
            array=np.full(n_obs,exptime,dtype=np.float64),
        ),
        fits.Column(
            name="ra",
            format="D",
            array=np.full(n_obs,9.55,dtype=np.float64),
        ),
        fits.Column(
            name="dec",
            format="D",
            array=np.full(n_obs,-44.1,dtype=np.float64),
        ),
        fits.Column(
            name="pa",
            format="D",
            array=np.full(n_obs,180,dtype=np.float64),
        ),
        fits.Column(
            name="filter",
            format="4A",
            array=np.full(n_obs,"H158"),
        ),
    ]
    
    table_hdu = fits.BinTableHDU.from_columns(columns)
    fits.HDUList(
        [
            fits.PrimaryHDU(),
            table_hdu,
        ]
    ).writeto(
        obsfile,
        overwrite=True)
    
    return obsfile

@pytest.fixture
def asdf_files(tmp_path):
    """
    """
    nside = 4088
    sca = 17
    ov = 1
    n = 121
    
    l2_dir = tmp_path/"L2"
    l21_dir = tmp_path/"L2_1"
    l2_dir.mkdir()
    l21_dir.mkdir()
    wcsobj = make_wcs()
    
    #compute psf once
    psf = PolychromaticPSF(sca,2043.5,2043.5,np.array([1.58]), frame="science").compute_poly_psf(
        postage_stamp_size=n, cycle=10, ovsamp=ov, use_filter="H")

    # desired star positions in (x,y) for GWCS
    detector_locations = {
        "star_1": (450.0,450.0),
        "star_2": (1050.0,950.0),
        "star_3": (250.0,750.0),
        "star_too_faint": (1234.0,1234.0),
        "star_too_old": (950.0,1550.0) }
    #convert (x,y) to (ra,dec)
    sky_positions = {}
    for star, (x,y) in detector_locations.items():
        ra_star,dec_star= wcsobj(x,y)
        sky_positions[star] = (float(ra_star),float(dec_star))

    #define fluxes for candidate stars as well as too faint ones
    bright_flux = 1.0e4
    faint_flix = 1.0e2

    #generate the sources that will appear in the asdf file, assuming starting at OBSID 9
    # with the given exptime, this will make obsid 2 too old to be considered
    stars_by_obsid = {
        #old source
        2: [("star_too_old",bright_flux)],
        3: [("star_1",bright_flux)],
        4: [("star_too_faint",faint_flux)],
        5: [("star_2",bright_flux)],
        6: [("star_3",bright_flux)],
        7: [("star_1",bright_flux)],
        8: [("star_3",bright_flux)],
    }
    #create all asdf files
    injected_pixels = {}
    for obsid in range(10):
        image = np.zeros((nside,nside),dtype=np.float32)
        injected_pixels[obsid] = {}
        for star, flux in stars_by_obsid.get(obsid,[]):
            ra_star,dec_star = sky_positions[star]
            ypix,xpix = inject_star(image=image,wcsobj=wcsobj,
                                    psf=psf,ra=ra_star,dec=dec_star,
                                    flux=flux,ov=ov,n=n)
            dq = np.zeros((nside,nside),dtype=np.uint32)
            dq |= np.where(image > 500.0, pixel.SATURATED, 0)
            es = np.where(image > 500.0, 6, -1)
            
            tree = {
                "processinfo": {
                    "endslice": es,
                },
                "roman": {
                    "meta": {
                        "wcs": wcsobj,
                        "instrument": {
                            "optical_element": "F158",
                            "detector": f"WFI{sca:02d}",
                        },
                    },
                    "data": image,
                    "dq": dq,
                },
            }
            outfile = ( l2_dir/f"sim_L2_H158_{obsid}_{sca}.asdf")
            asdf.AsdfFile(tree).write_to(outfile,all_array_compression="zlib")
        # generate one asdf file in the L2.1 directory to get the ball rolling
        current_image = np.zeros((nside,nside),dtype=np.float32)
        current_dq = np.zeros((nside,nside),dtype=np.uint32)
        current_es = np.full((nside,nside),-1,dtype=np.int8)
        current_tree = {
            "processinfo": {
                "endslice": current_es,
            },
            "roman": {
                "meta": {
                    "wcs": wcsobj,
                    "instrument": {
                        "optical_element": "F158",
                        "detector": f"WFI{sca:02d}",
                    },
                },
                "data": current_image,
                "dq": current_dq,
            },
        }
        current_path = (l21_dir/f"sim_L2_1_H158_9_{sca}.asdf")
        asdf.AsdfFile(current_tree).write_to(current_path,all_array_compression="zlib")
        
        return {
            "l2_dir": l2_dir,
            "l21_dir": l21_dir,
            "current_path": current_path,
            "pixels": {
                name: (int(round(y)),int(round(x)))
                for name, (x,y) in detector_positions.items()
            },
            "expected_obsids": [8,7,6,5,4,3]
        }


@pytest.fixture
def config(tmp_path,obsfile,asdf_files):
    """
    """
    l21_dir = asdf_files["l21_dir"]
    l21_path = f"{l21_dir}/"
    config = {
        "OBSFILE": str(obsfile),
        "INDATA": [l21_path,"L2_2506"],
        "TILESCHM": "FourSquare-Dec2025",
        "RERUN": "11",
        "MOSAIC": 1,
        "CTR": [
    	9.55,
    	-44.1
        ],
        "LONPOLE": 180.0,
        "OUTSIZE": [
            60,
    	34,
    	0.049019607843137254
        ],
        "BLOCK": 40,
        "FILTER": 2,
        "LAKERNEL": "Cholesky",
        "KAPPAC": [
             5e-4
        ],
        "INPSF": [
    	"/fs/scratch/PCON0003/cond0007/Jul26/psf",
            "L2_2506",
            6
        ],
        "EXTRAINPUT": [
            "truth,0.004906087669824225",
            "gsstar14",
            "nstar14,2e5,86,3",
            "gsext14,seed=100",
            "gsext14,seed=100,shear=0.02:0",
            "gsext14,seed=100,shear=0:0.02",
            "gsext14,seed=100,shear=-0.02:0",
            "gsext14,seed=100,shear=0:-0.02",
            "gsext14,seed=100,shape=0:0",
            "gsext14,seed=100,shape=0:0,shear=0.02:0",
            "gsext14,seed=100,shape=0:0,shear=0:0.02",
            "gsext14,seed=100,shape=0:0,shear=-0.02:0",
            "gsext14,seed=100,shape=0:0,shear=0:-0.02",
            "gsext14,seed=100,hlr=0.2,shape=0:0",
            "gsext14,seed=100,hlr=0.2,shape=0:0,shear=0.02:0",
            "gsext14,seed=100,hlr=0.2,shape=0:0,shear=0:0.02",
            "gsext14,seed=100,hlr=0.2,shape=0:0,shear=-0.02:0",
            "gsext14,seed=100,hlr=0.2,shape=0:0,shear=0:-0.02",
            "whitenoise10",
            "noise,Rz4PbrS2C1",
            "noise,Rz4PbrS2C2",
            "noise,Rz4PbrS2C3",
            "noise,Rz4PbrS2C4",
            "noise,Rz4OS2C5",
            "noise,Rz4OS2C6",
            "noise,Rz4OS2C7",
            "noise,Rz4OS2C8"
        ],
        "PADSIDES": "all",
        "OUTMAPS": "USTN",
        "OUT": "/fs/scratch/PCON0003/cond0007/Jul26/H1/coadd/H1-r11",
        "INPAD": 0.70,
        "NPIXPSF": 18,
        "FADE": 1,
        "PAD": 1,
        "NOUT": 1,
        "OUTPSF": "GAUSSIAN",
        "EXTRASMOOTH": 0.9265328730414752,
        "TEMPFILE": "/fs/scratch/PCON0003/cond0007/tmp/pyimcomrun_2X",
        "INLAYERCACHE": "/fs/scratch/PCON0003/cond0007/temp/cache2/r2",
        "PSFINTERP": "G4460",
        "PSFSPLIT": [5.25, 8.75, 0.01, true],
        "DSMODEL": ["constant",
            4088],
        "DSOBSFILE": "/fs/scratch/PCON0003/cond0007/temp/fits-H/tempimage_",
        "DSOUT": ["/fs/scratch/PCON0003/cond0007/temp/ds-H/",
            "_ds_out.txt"
        ],
        "CGMODEL": ["PR",
            3,
            1e-3
            ],
        "DSCOST": ["quadratic",
            0,
            0
            ]
    }
    config_path = (tmp_path/"persistence_test_config.json")
    with config_path.open("w",encoding="utf-8") as config_file:
        json.dump(config,config_file,indent=4)
        
    return config_path


@pytest.fixture
def observation_table(tmp_path):
    """
    Create a temporary FITS table whose observations are deliberately
    not stored in chronological order.

    Sort Order   Python index    DATE
    ----------   ------------    -----------
       2              0          62000.03000
       0              1          62000.01000
       3              2          62000.04000
       1              3          62000.02000
    Chronological order:
        row 1 -> row 3 -> row 0 -> row 2
    """
    dates = np.array([62000.03000, 62000.01000, 62000.04000, 62000.02000], dtype=np.float64)

    date_column = fits.Column(
        name="date",
        format="D",
        array=dates,
    )

    primary_hdu = fits.PrimaryHDU()
    table_hdu = fits.BinTableHDU.from_columns([date_column])

    filepath = tmp_path / "observation_table.fits"
    fits.HDUList([primary_hdu, table_hdu]).writeto(filepath)

    return filepath


def test_previous_obsid_exact_match(observation_table):
    """An exact OBSID match should return the next-lowest OBSID."""
    result = previous_obsid(
        str(observation_table),
        obsid=62000.03000,
    )

    assert result == pytest.approx(62000.02000)


def test_previous_obsid_between_observations(observation_table):
    """
    An OBSID between two recorded observations should return the
    largest recorded OBSID below it.
    """
    result = previous_obsid(
        str(observation_table),
        obsid=62000.03500,
    )

    assert result == pytest.approx(62000.03000)


def test_previous_obsid_after_last_observation(observation_table):
    """An OBSID after every observation should return the final OBSID."""
    result = previous_obsid(
        str(observation_table),
        obsid=62000.05000,
    )

    assert result == pytest.approx(62000.04000)


@pytest.mark.parametrize(
    "obsid",
    [
        62000.01000,  # Exactly the earliest observation
        62000.00000,  # Earlier than every observation
    ],
)
def test_previous_obsid_returns_none_when_no_previous_exists(
    observation_table,
    obsid,
):
    """Checks that None is passed if the observation_table does
    not contain an obsid previous to the one passed.
    """
    result = previous_obsid(
        str(observation_table),
        obsid=obsid,
    )

    assert result is None


def test_previous_row_number(observation_table):
    """
    As passed, row_number = 3 corresponds to DATE=62000.02000.

    The previous chronological observation is DATE=62000.01000,
    which is located in FITS row 1.
    """
    result = previous_obsid(
        str(observation_table),
        row_number=3,
    )

    assert result == 1


def test_first_chronological_row_has_no_previous_row(observation_table):
    """
    FITS row 1 contains the earliest observation, so it has no
    chronologically previous observation.
    """
    result = previous_obsid(
        str(observation_table),
        row_number=1,
    )

    assert result is None


def test_passing_both_search_arguments_prints_warning(
    observation_table,
):
    """Checks that ValueError is raised if we provide both arguments."""
    with pytest.raises(ValueError):
        previous_obsid(
            str(observation_table),
            row_number=2,
            obsid=62000.03000,
        )


def test_passing_no_search_argument_returns_none(observation_table):
    """Checks to make sure None is passed when neither
    row_number or obsid is passed
    """
    result = previous_obsid(str(observation_table))

    assert result is None


def test_get_prev_obs_returns_none_when_no_previous_exists(observation_table):
    """Checks that Nones are passed if the observation_table does
    not contain an id previous to the one passed.
    """
    result1, result2 = get_prev_obs(str(observation_table),id=1)

    assert result1 is None and result2 is None


def test_get_prev_obs_nbackup_not_supplied(observation_table):
    """
    As passed, row_number = 2 corresponds to DATE=62000.04000.

    The previous chronological observation is DATE=62000.03000,
    which is located in FITS row 0, and their difference is 0.01.
    """
    result1, result2 = get_prev_obs(str(observation_table), id=2)

    assert result1 == 0 and np.abs(result2 - 0.01) < 1.0e-5


def test_get_prev_obs_with_nbackup_passed(observation_table):
    """
    As passed, row_number = 0 corresponds to DATE=62000.03000.

    The observation that is (n=2) before is DATE=62000.01000,
    which is located in FITS row 1, and their difference is 0.02.
    """
    result1, result2 = get_prev_obs(str(observation_table), id=0, nbackup=2)

    assert result1 == 1 and np.abs(result2 - 0.02) < 1.0e-5


#make a pytest.fixture here for: writing a config json, 
# make a sample observation table fits,
# (can use OU24 tables for formatting), and make a
# series of corresponding asdf images (can use the template of test_many_stars)
# and re-use the same PSF over and over instead of re-drawing it for each sca
# then,
def test_persistence_flagging(config,asdf_files):
    """ Make sure persistence.run() returns a persistence mask and
    a list of the correct obsids used to build it, namely: (INSERT
    APPLICABLE OBSIDS)
    """
    mask_list = run(str(config))
    
    assert isinstnace(mask_list,list)
    assert len(mask_list) == 1
    
    persistence_mask, prev_obsids = mask_list[0]
    
    assert isinstance(persistence_mask,np.ndarray)
    assert persistence_mask.shape == (4088,4088)
    assert prev_obsids == [8,7,6,5,4,3]
    
    pixels = asdf_files["pixels"]
    
    # all of these should be considered
    assert persistence_mask[pixels["star_1"]]
    assert persistence_mask[pixels["star_2"]]
    assert persistence_mask[pixels["star_3"]]
    
    #this one should be too faint
    assert not persistence_mask[pixels["star_too_faint"]]
    
    #this one should be too old
    assert not persistence_mask[pixels["star_too_bright"]]
