"""
Test functions for internal object finding.
"""

import asdf
import numpy as np
from astropy import coordinates, units
from astropy.modeling import models
from gwcs import coordinate_frames, wcs
from psfsim.polychrom import PolychromaticPSF
from roman_datamodels.dqflags import pixel
from roman_hlis_l2_driver.starcatalogs.internal import brightobj_from_manyimg, encirc_center
from roman_hlis_l2_driver.starcatalogs.interp import extract_spikedata, interp5pt
from roman_hlis_l2_driver.starcatalogs.main import spikedata_to_mask, stardata_to_mask_wcs
from scipy.ndimage import shift


def test_spikedata_to_mask():
    """Test for converting spike information into a 2D mask."""

    sp_theta = np.linspace(15, 345, 12) * np.pi / 180.0
    sp_length = np.array([60, 60, 60, 60, 60, 100, 60, 60, 60, 60, 60, 60])
    sp_width = np.array([6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 3])

    mask1 = spikedata_to_mask(4088, 2000.0, 45.0, sp_theta, sp_length, sp_width, np.pi / 180.0 * 6.5)

    assert np.all(~mask1[105:, :])
    assert np.all(mask1[104, 1986:1992] == np.array([False, True, True, True, True, False]))
    assert np.all(mask1[104, 2009:2015] == np.array([False, True, True, True, True, False]))
    assert np.count_nonzero(mask1[:, 1922]) == 13
    assert np.count_nonzero(mask1[25:45, 2039]) == 8


def test_interp5pt():
    """Interpolation test function."""

    u = 0.7
    v = 0.2
    upos = np.array([-1, 1, 0, -1, 1, u])
    vpos = np.array([-1, -1, 0, 1, 1, v])
    image = np.zeros((40, 6))
    for i in range(40):
        image[i, :] = (
            np.cos(i) * upos
            - i * np.sin(i) * vpos
            + 0.9**i * upos * vpos
            + np.exp(-np.sin(2 * i))
            + (i // 10) * upos * vpos
        )
    i2 = image[:, :-1]
    assert np.allclose(image[:, -1], interp5pt(i2, u, v, 1))

    # now try a few other orientations
    assert np.allclose(image[:, -1], interp5pt(i2.T, u, v, 0))
    assert np.allclose(image[:, -1].reshape((5, 8)), interp5pt(i2.reshape((5, 8, 5)), u, v, 2))

    # just checking this isn't all 0's so it didn't fail trivially
    assert np.all(np.abs(i2) >= 0.01)


def test_extract_spikedata():
    """Simple test for getting diffraction spike information."""

    ls = [14, 22, 64, 100]
    ws = [2.2, 2.6, 4.4, 5.5]
    for j in range(4):
        d = extract_spikedata("H", 1, 4087.5, 4087.5, 3e-4 / 10**j)
        print(d)
        assert np.all(np.abs(np.degrees(d[0]) - np.linspace(15, 345, 12)) < 10)
        assert np.all(np.abs(np.log(d[1] / ls[j])) < 0.2)
        assert np.all(np.abs(np.log(d[2] / ws[j])) < 0.2)


def test_encirc_center():
    """Test function for encirc_center."""

    x_, y_ = np.meshgrid(np.arange(7), np.arange(8))
    use = (x_ - 3) ** 2 + (y_ - 3) ** 2 < 10.5
    use[3, 3:5] = False
    use[6, 5] = True
    use[0, 2] = False
    use[0, 4] = False
    print(use)

    # 0 values
    use2 = np.zeros_like(use)
    x, y, r = encirc_center(use2)
    assert abs(x + 1) < 1e-6
    assert abs(y + 1) < 1e-6
    assert abs(r + 1) < 1e-6
    print(x, y, r)

    # 1 value
    use2[5, 2] = True
    x, y, r = encirc_center(use2)
    assert abs(x - 2) < 1e-6
    assert abs(y - 5) < 1e-6
    assert abs(r) < 1e-6
    print(x, y, r)

    # 2 values
    use2[3, 3] = True
    x, y, r = encirc_center(use2)
    assert abs(x - 2.5) < 1e-6
    assert abs(y - 4) < 1e-6
    assert abs(r - 1.118033988749895) < 1e-6
    print(x, y, r)

    # General
    x, y, r = encirc_center(use)
    print(x, y, r)
    assert abs(x - 3.0454545454545454) < 1e-6
    assert abs(y - 3.3181818181818183) < 1e-6
    assert abs(r - 3.3184931360807237) < 1e-6


def test_mask_many(tmp_path):
    """Test star mask functions for many images."""

    tmp_path = str(tmp_path)
    # list of obsid, sca, ra, dec, lonpole for this test
    obslist = [
        (145, 11, 16.5, -19.5, 180.0),
        (148, 12, 16.5, -19.4, 225.0),
    ]

    # starlist: ra, dec, flux
    stars = [
        (16.5, -19.36, 100.0),
        (16.515, -19.52, 1.0e5),
        (16.49, -19.438, 2.0e5),
        (16.495, -19.438, 4.0e5),
    ]

    # helpful
    deg = np.pi / 180.0
    nside = 4088
    nbdy = 256  # boundary pixels
    ov = 4  # oversampling factor for PSF
    n = 121  # size of postage stamp to draw

    # Loop over observations
    for obs in obslist:
        (obsid, sca, ra, dec, lonpole) = obs

        # Construct the GWCS
        shiftscale = models.Shift(-2043.5) | models.Scale(0.11 / 3600.0)
        det2sky = (
            ((shiftscale | models.Scale(-1)) & shiftscale)
            | models.Pix2Sky_TAN()
            | models.RotateNative2Celestial(ra, dec, lonpole)
        )
        wcsobj = wcs.WCS(
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

        # include boundary, trim later
        extended_image = np.zeros((nside + 2 * nbdy, nside + 2 * nbdy), dtype=np.float32)

        for s in stars:
            (ra_s, dec_s, flux_s) = s

            # weed out stars more than 0.1 deg from the SCA center
            mu = np.sin(dec_s * deg) * np.sin(dec * deg) + np.cos(dec_s * deg) * np.cos(dec * deg) * np.cos(
                (ra_s - ra) * deg
            )
            if mu < np.cos(0.1 * deg):
                continue

            # where is this?
            x, y = wcsobj.invert(ra_s, dec_s)
            xc = int(np.round(x))
            yc = int(np.round(y))
            print(x, y, ra_s, dec_s, flux_s)
            nbdy2 = nbdy - n // 2
            if xc < -nbdy2 or xc >= nside + nbdy2 or yc < -nbdy2 or yc >= nside + nbdy2:
                continue  # skip if not in boundary

            psf = PolychromaticPSF(sca, x, y, np.array([1.58]), frame="science").compute_poly_psf(
                postage_stamp_size=n, cycle=10, ovsamp=ov, use_filter="H"
            )
            stamp = np.sum(
                shift(psf, (ov * (x - xc), ov * (y - yc)), order=1, mode="constant").reshape((n, ov, n, ov)),
                axis=(1, 3),
            )
            extended_image[
                yc + nbdy - n // 2 : yc + nbdy + n // 2 + 1, xc + nbdy - n // 2 : xc + nbdy + n // 2 + 1
            ] += flux_s * stamp

        # trim extension region
        image = extended_image[nbdy:-nbdy, nbdy:-nbdy]

        # set data quality flags
        dq = np.zeros((nside, nside), dtype=np.uint32)
        dq |= np.where(image > 500.0, pixel.SATURATED, 0)
        es = np.where(image > 500.0, 6, -1)

        # save as ASDF
        asdf.AsdfFile(
            {
                "processinfo": {"endslice": es},
                "roman": {
                    "meta": {
                        "wcs": wcsobj,
                        "instrument": {"optical_element": "F158", "detector": f"WFI{sca:02d}"},
                    },
                    "data": image,
                    "dq": dq,
                },
                "mask": np.zeros((4088, 4088), dtype=np.uint8),
            }
        ).write_to(f"{tmp_path}/obs_{obsid}_{sca}.asdf")

        from astropy.io import fits

        fits.PrimaryHDU(image).writeto(f"{tmp_path}/obs_{obsid}_{sca}.fits")

    # Now run the star finder
    stars_recovered = brightobj_from_manyimg(tmp_path + "/obs_{0:d}_{1:d}.asdf")
    print(stars_recovered)
    assert len(stars_recovered) == 5
    print("\n")

    stars_recovered = brightobj_from_manyimg(tmp_path + "/obs_{0:d}_{1:d}.asdf", clean=True, verbose=True)
    stars_unique = stars_recovered[stars_recovered["nchild"] > 0]
    assert np.all(stars_unique["idsca"] == np.array([14812, 14812]))
    assert np.all(np.abs(stars_unique["ra"] - np.array([16.49, 16.495])) < 5.0e-5)
    assert np.all(np.abs(stars_unique["dec"] - np.array([-19.438, -19.438])) < 5.0e-5)
    assert np.all(np.abs(np.log(stars_unique["eflux"] / np.array([2.0e5, 4.0e5]))) < 0.25)

    n_gt2 = [3032, 2483]

    for obs in obslist:
        (obsid, sca, ra, dec, lonpole) = obs
        l2 = f"{tmp_path}/obs_{obsid}_{sca}.asdf"
        thismask = stardata_to_mask_wcs(l2, stars_unique, 0.1)
        if (obsid, sca) == (145, 11):
            assert np.all(~thismask[:2044, :])
            assert np.count_nonzero(thismask[3964, 2100:2200]) == 14
            assert np.count_nonzero(thismask[3964, 2200:2250]) == 11
            assert np.count_nonzero(thismask[4030:4070, 2287]) == 16
            assert np.all(~thismask[4069:4076, 2349:2356])

        with asdf.open(l2) as a:
            im = np.where(~thismask, a["roman"]["data"], -1)
        assert np.abs(np.count_nonzero(im > 2) - n_gt2[0]) < 10

        # uncomment if you want to make images for debugging
        # from astropy.io import fits
        # fits.PrimaryHDU(im).writeto(f"{tmp_path}/obs_{obsid}_{sca}_masked.fits")

        # remove first item
        n_gt2 = n_gt2[1:]
