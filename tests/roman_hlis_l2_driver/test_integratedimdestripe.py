"""Integrated test for imdestripe."""

import gzip
import os
import shutil
from urllib.request import urlretrieve

import asdf
import numpy as np
import pytest
from astropy.io import fits
from roman_hlis_l2_driver.destripe_interface import destripe, destripe_setup
from roman_hlis_l2_driver.outliers import outlier_flagging

# Example configuration file.
# Note that we will replace $TMPDIR.
# The only subfolders that are actually used in this test are:
#
# - $TMPDIR/L2/
# - $TMPDIR/fits-F/
# - $TMPDIR/ds-F/

EXAMPLE_CONFIG = """{
    "OBSFILE": "$TMPDIR/Roman_WAS_obseq_11_1_23.fits",
    "INDATA": [
        "$TMPDIR/L2/",
        "L2_2506"
    ],
    "TILESCHM": "FourSquare-Dec2025",
    "RERUN": "06",
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
    "FILTER": 1,
    "LAKERNEL": "Cholesky",
    "KAPPAC": [
         5e-4
    ],
    "INPSF": [
	"$TMPDIR/psf",
        "L2_2506",
        6
    ],
    "EXTRAINPUT": [
        "truth,0.004906087669824225",
        "gsstar14"
    ],
    "PADSIDES": "all",
    "OUTMAPS": "USTN",
    "OUT": "$TMPDIR/imtest-F1",
    "INPAD": 0.70,
    "NPIXPSF": 18,
    "FADE": 1,
    "PAD": 1,
    "NOUT": 1,
    "OUTPSF": "GAUSSIAN",
    "EXTRASMOOTH": 0.9265328730414752,
    "TEMPFILE": "$TMPDIR/pyimcomrun_2X",
    "INLAYERCACHE": "$TMPDIR/r4",
    "PSFINTERP": "G4460",
    "PSFSPLIT": [5.25, 8.75, 0.01],
    "DSMODEL": ["constant",
        4088],
    "DSOBSFILE": "$TMPDIR/fits-F/tempimage_",
    "DSOUT": ["$TMPDIR/ds-F/",
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
"""


def test_null(tmp_path):
    """Tests when you give destripe_one_layer zero files."""

    tmp_path = str(tmp_path)  # convert to string

    # first, get the configuration file.
    with open(tmp_path + "/cfg.txt", "w") as f:
        f.write(EXAMPLE_CONFIG.replace("$TMPDIR", str(tmp_path)))

    # now make a directory
    os.makedirs(tmp_path + "/L2", exist_ok=True)
    assert destripe.destripe_one_layer(tmp_path + "/cfg.txt") is None


def _collect(files, newloc):
    """Collects files for the test, sends to `newloc`."""

    for f in files:
        urlretrieve(
            "https://github.com/Roman-HLIS-Cosmology-PIT/pyimcom/wiki/test-files/imdestripe/" + f,
            newloc + "/" + f,
        )
        if f[-3:] == ".gz":
            f_final = newloc + "/" + f[:-3]
            with gzip.open(f_final + ".gz", "rb") as f1, open(f_final, "wb") as f2:
                shutil.copyfileobj(f1, f2)
            os.remove(f_final + ".gz")


def test_integrated(tmp_path):
    """Integrated test for imdestripe."""

    tmp_path = str(tmp_path)  # convert to string

    # first, get the configuration file.
    with open(tmp_path + "/cfg.txt", "w") as f:
        f.write(EXAMPLE_CONFIG.replace("$TMPDIR", str(tmp_path)))

    # now download the files that imdestripe needs
    os.makedirs(tmp_path + "/L2", exist_ok=True)
    cpath = [
        "sim_L2_F184_1433_11.asdf",
        "sim_L2_F184_14844_6.asdf",
        "sim_L2_F184_1433_12.asdf",
        "sim_L2_F184_1433_11_noise.asdf",
        "sim_L2_F184_14844_6_noise.asdf",
        "sim_L2_F184_1433_12_noise.asdf",
        "sim_L2_F184_1433_11_mask.fits.gz",
        "sim_L2_F184_14844_6_mask.fits.gz",
        "sim_L2_F184_1433_12_mask.fits.gz",
    ]
    _collect(
        cpath,
        tmp_path + "/L2",
    )

    # test moving files as if they are noise
    os.makedirs(tmp_path + "/asifnoise-F", exist_ok=True)
    destripe_setup.setup_all_files(
        tmp_path + "/L2/sim_L2_F184", tmp_path + "/asifnoise-F/nf", noiseid=0, verbose=True
    )
    with asdf.open(tmp_path + "/L2/sim_L2_F184_1433_11.asdf") as a1:
        block = np.copy(a1["roman"]["data"][:16, :16])
    with asdf.open(tmp_path + "/L2/sim_L2_F184_1433_11_noise.asdf") as a2:
        block += a2["noise"][0, :16, :16]
    with fits.open(tmp_path + "/asifnoise-F/nf_1433_11.fits") as f:
        assert np.allclose(f[0].data[:16, :16], block, atol=1e-5, rtol=1e-5)
    # now clear old files (this part also asserts that they exist!)
    delfiles = [
        "asifnoise-F/nf_14844_6.fits",
        "asifnoise-F/nf_1433_11.fits",
        "asifnoise-F/nf_1433_12.fits",
    ]
    for df in delfiles:
        os.remove(tmp_path + "/" + df)

    # make directories for imdestripe
    os.makedirs(tmp_path + "/ds-F", exist_ok=True)
    os.makedirs(tmp_path + "/fits-F", exist_ok=True)

    # get parameters before destriping
    obsid = [1433, 1433, 14844]
    sca = [11, 12, 6]
    iqr_old = []
    for j in range(3):
        with asdf.open(tmp_path + f"/L2/sim_L2_F184_{obsid[j]:d}_{sca[j]:d}.asdf") as f:
            rows = np.median(f["roman"]["data"], axis=1)
        iqr_old.append(np.percentile(rows, 75) - np.percentile(rows, 25))
    for i in range(2):
        for j in range(3):
            with asdf.open(tmp_path + f"/L2/sim_L2_F184_{obsid[j]:d}_{sca[j]:d}_noise.asdf") as f:
                rows = np.median(f["noise"][i, :, :], axis=1)
            iqr_old.append(np.percentile(rows, 75) - np.percentile(rows, 25))
    iqr_old = np.array(iqr_old)
    print("IQR OLD =", iqr_old)

    # run the destriping
    destripe.destripe_all_layers(tmp_path + "/cfg.txt", verbose=True, tracktest=True)

    # check overlap matrix
    mtarget = np.array([[1.0, 0.0, 0.05942073], [0.0, 1.0, 0.33928332], [0.05942073, 0.33928332, 1.0]])
    mat = np.load(tmp_path + "/ds-F/ovmat.npy")
    assert np.all(np.abs(mat - mtarget) < 0.02)

    # check residuals
    resids = np.loadtxt(tmp_path + "/ds-F/cg_log.csv", delimiter=",", skiprows=1)[:, 1]
    targets = np.array([4.34582306e03, 1.12769622e03, 4.80621360e02])
    assert np.all(resids / targets > 0.8)
    assert np.all(resids / targets < 1.2)

    # checks of parameters
    with fits.open(tmp_path + "/ds-F/final_params.fits") as f:
        pars = f[0].data
    print(np.std(pars, axis=1))
    assert np.std(pars[0, :]) < 0.002
    assert 0.006 < np.std(pars[1, :]) < 0.01
    assert 0.006 < np.std(pars[2, :]) < 0.01

    iqr_new = []
    for j in range(3):
        with asdf.open(tmp_path + f"/L2/sim_L2_F184_{obsid[j]:d}_{sca[j]:d}.asdf") as f:
            rows = np.median(f["roman"]["data"], axis=1)
        iqr_new.append(np.percentile(rows, 75) - np.percentile(rows, 25))
    for i in range(2):
        for j in range(3):
            with asdf.open(tmp_path + f"/L2/sim_L2_F184_{obsid[j]:d}_{sca[j]:d}_noise.asdf") as f:
                rows = np.median(f["noise"][i, :, :], axis=1)
            iqr_new.append(np.percentile(rows, 75) - np.percentile(rows, 25))
    iqr_new = np.array(iqr_new)
    print("IQR NEW =", iqr_new)
    ratio = iqr_new / iqr_old
    print("RATIO =", ratio)
    ratio_target = np.array(
        [1.0, 0.49180344, 0.46751416, 1.0, 0.48641852, 0.42708334, 1.0, 0.44493178, 0.4264151]
    )
    assert np.all(np.abs(ratio - ratio_target) < 0.2)

    # now run outliers!
    # this test won't actually flag anything because of insufficient overlap, but it exercises the mechanics
    om = outlier_flagging.OutlierMap(
        tmp_path + "/cfg.txt",
        max_workers=2,
        run_and_save=tmp_path + "/mask_F1.asdf",
        mask_kwargs={"cut_c": 0.3},
    )
    with asdf.open(tmp_path + "/mask_F1.asdf") as a:
        assert a["N"] == 3
        assert len(a["files"]) == 3
        assert a["files"][0][-24:] == "sim_L2_F184_1433_11.asdf"
        assert np.shape(a["mask"]) == (3, 4088, 4088)
        assert np.all(np.abs(a["overlap"] - mtarget) < 0.02)
    # exercise some other OutlierMap functionality
    with pytest.raises(ValueError):
        om.save_data("haha.nosuchextension")
    assert om._worker_outlier_mask(0, {"cut_c": 0.3})[0] == 0
    om2 = outlier_flagging.OutlierMap(tmp_path + "/cfg.txt", max_files=1)
    assert om2.n_obs == 1

    # now clear old files (this part also asserts that they exist!)
    delfiles = [
        "ds-F/F184_DS_14844_6.fits",
        "ds-F/F184_DS_1433_12.fits",
        "ds-F/F184_DS_1433_11.fits",
        "ds-F/final_params.fits",
        "ds-F/cg_log.csv",
        "ds-F/ovmat.npy",
        "fits-F/tempimage_F184_1433_12.fits",
        "fits-F/tempimage_F184_14844_6.fits",
        "fits-F/tempimage_F184_1433_11.fits",
    ]
    for df in delfiles:
        os.remove(tmp_path + "/" + df)
    for df in cpath:
        df2 = df[:-3] if df[-3:] == ".gz" else df
        os.remove(tmp_path + "/L2/" + df2)


if __name__ == "__main__":
    test_integrated("out")
