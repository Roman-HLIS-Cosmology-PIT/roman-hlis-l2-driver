"""Test functions for full field of view images."""

import os

import asdf
import gwcs
import numpy as np
import pytest
from astropy import coordinates as coord
from astropy import units as u
from astropy import wcs
from astropy.io import fits
from astropy.modeling import models
from gwcs import coordinate_frames as cf
from numpy.random import RandomState
from roman_hlis_l2_driver.fullfovexport.fullfov import FullFoVImage, FullFoVImageFromFile

# CD and CRPIX values for linear approximations to the SCA WCSs.
wcsdata = np.array(
    [
        [
            3.08219178082192e-05,
            1.22309197651664e-07,
            0,
            -3.02103718199609e-05,
            9.31141597190574e-10,
            -254.695456590193,
            819.255060728745,
        ],
        [
            3.04549902152642e-05,
            1.22309197651664e-07,
            -1.22309197651663e-07,
            -2.97211350293542e-05,
            9.05141916965698e-10,
            -302.076620500445,
            5721.07850461111,
        ],
        [
            3.02103718199609e-05,
            2.44618395303327e-07,
            -1.22309197651663e-07,
            -2.88649706457926e-05,
            8.71991576701989e-10,
            -340.491336421342,
            10368.6800480357,
        ],
        [
            3.05772994129159e-05,
            1.22309197651664e-07,
            -2.44618395303327e-07,
            -3.02103718199609e-05,
            9.23721665434799e-10,
            -4684.7680248753,
            -19.9937811750988,
        ],
        [
            3.04549902152642e-05,
            3.6692759295499e-07,
            -3.66927592954991e-07,
            -2.97211350293542e-05,
            9.05022240647057e-10,
            -4754.73767727858,
            4920.56054745611,
        ],
        [
            3.02103718199609e-05,
            4.89236790606654e-07,
            -3.66927592954991e-07,
            -2.91095890410959e-05,
            8.79231993979803e-10,
            -4894.76339878177,
            9448.63987477456,
        ],
        [
            3.04549902152642e-05,
            2.44618395303327e-07,
            -3.66927592954991e-07,
            -3.03326810176125e-05,
            9.23691746355138e-10,
            -9119.77676286723,
            -2073.79302302983,
        ],
        [
            3.02103718199609e-05,
            4.89236790606654e-07,
            -3.66927592954991e-07,
            -2.98434442270059e-05,
            9.0140203200815e-10,
            -9255.53159851301,
            2786.07620817844,
        ],
        [
            2.99657534246575e-05,
            6.11545988258318e-07,
            -6.11545988258311e-07,
            -2.92318982387476e-05,
            8.75581866261235e-10,
            -9476.57488467452,
            7313.76934905176,
        ],
        [
            3.06996086105675e-05,
            -1.22309197651664e-07,
            1.22309197651664e-07,
            -3.02103718199609e-05,
            9.27431631312686e-10,
            4156.44544809343,
            827.807471449771,
        ],
        [
            3.04549902152642e-05,
            -2.44618395303327e-07,
            1.22309197651664e-07,
            -2.95988258317025e-05,
            9.0140203200815e-10,
            4207.94795539033,
            5735.52044609665,
        ],
        [
            3.02103718199609e-05,
            -2.44618395303327e-07,
            1.22309197651664e-07,
            -2.88649706457926e-05,
            8.71991576701989e-10,
            4262.97958483445,
            10367.9787270544,
        ],
        [
            3.05772994129158e-05,
            -2.44618395303327e-07,
            2.44618395303327e-07,
            -3.00880626223092e-05,
            9.1995186139759e-10,
            8568.20762326005,
            -30.0470924938209,
        ],
        [
            3.03326810176125e-05,
            -3.6692759295499e-07,
            2.44618395303327e-07,
            -2.97211350293542e-05,
            9.0143195108781e-10,
            8671.99004281589,
            4891.17687278038,
        ],
        [
            3.02103718199609e-05,
            -3.66927592954991e-07,
            3.66927592954991e-07,
            -2.89872798434442e-05,
            8.75581866261235e-10,
            8754.5221937468,
            9476.99395181958,
        ],
        [
            3.03326810176125e-05,
            -2.44618395303327e-07,
            2.44618395303327e-07,
            -3.03326810176125e-05,
            9.2001169955691e-10,
            13087.5824390244,
            -2119.77756097561,
        ],
        [
            3.03326810176125e-05,
            -3.66927592954991e-07,
            4.89236790606654e-07,
            -2.99657534246575e-05,
            9.08762125604605e-10,
            13130.6172384276,
            2825.69171001514,
        ],
        [
            3.00880626223092e-05,
            -6.11545988258318e-07,
            4.89236790606654e-07,
            -2.92318982387476e-05,
            8.79231993979802e-10,
            13350.5118589853,
            7261.98346207507,
        ],
    ]
)


def make_simple_wcs(ra, dec, pa, sca):
    """
    Makes a simple approximate WCS for an SCA.

    Parameters
    ----------
    ra : float
        RA of the WFI center in degrees.
    dec : float
        Dec of the WFI center in degrees.
    pa : float
        PA of the WFI center in degrees.
    sca : int
        SCA number, in 1 .. 18.

    Returns
    -------
    astropy.wcs.WCS
        Simple WCS approximation.

    """

    outwcs = wcs.WCS(naxis=2)
    outwcs.wcs.crpix = [wcsdata[sca - 1, -2], wcsdata[sca - 1, -1]]
    outwcs.wcs.cd = wcsdata[sca - 1, :4].reshape((2, 2))
    outwcs.wcs.ctype = ["RA---ARC", "DEC--ARC"]
    outwcs.wcs.crval = [ra, dec]
    outwcs.wcs.lonpole = pa - 180.0 if pa >= 180.0 else pa + 180.0

    return outwcs


def test_ffov(tmp_path):
    """Test for the full field of view functions."""

    tmp_path = str(tmp_path)

    # this us just a dummy read pattern
    read_pattern = [[0], [1], [2, 3], [4, 5, 6, 7], [8, 9], [10]]
    frame_time = 3.15
    ngrp = len(read_pattern)
    tbar = np.zeros(ngrp, dtype=np.float32)
    tau = np.zeros(ngrp, dtype=np.float32)
    for i in range(ngrp):
        t0 = read_pattern[i][0]
        n = len(read_pattern[i])
        tbar[i] = frame_time * (t0 + (n - 1) / 2.0)
        tau[i] = frame_time * (t0 + (n - 1) * (2 * n - 1) / (6 * n))

    rs = RandomState(2026)  # legacy RNG

    # pointing
    ra = 150.0
    dec = -20.0
    pa = 60.0

    # first, make some dummy L2 files.
    for sca in range(1, 19):
        print("Making SCA", sca)

        # for SCA 11, let's make this file missing and see if the code can handle it
        if sca == 11:
            continue

        # make noise image
        image = rs.normal(loc=0.0, scale=0.35, size=(4088, 4088)).astype(np.float32)

        # now draw a line on it
        image[200:205, 400:800] += 1.0

        # now draw the SCA number in binary
        for j in range(5):
            b = 1 if (sca >> j) & 1 else -1
            image[100:120, 500 - 25 * j : 520 - 25 * j] += 10 * b

        # dq flags
        dq = np.zeros((4088, 4088), dtype=np.uint32)

        # pipeline version of the WCS
        this_wcs = make_simple_wcs(ra, dec, pa, sca)
        distortion = models.AffineTransformation2D(this_wcs.wcs.cd, translation=[0, 0])
        distortion.inverse = models.AffineTransformation2D(np.linalg.inv(this_wcs.wcs.cd), translation=[0, 0])
        celestial_rotation = models.RotateNative2Celestial(
            this_wcs.wcs.crval[0], this_wcs.wcs.crval[1], this_wcs.wcs.lonpole
        )
        det2sky = (
            (models.Shift(-(this_wcs.wcs.crpix[0] - 1)) & models.Shift(-(this_wcs.wcs.crpix[1] - 1)))
            | distortion
            | models.Pix2Sky_ARC()
            | celestial_rotation
        )
        det2sky.name = "mickey_mouse_invented_some_coordinate_system_haha"
        detector_frame = cf.Frame2D(name="detector", axes_names=("x", "y"), unit=(u.pix, u.pix))
        sky_frame = cf.CelestialFrame(reference_frame=coord.ICRS(), name="icrs", unit=(u.deg, u.deg))
        sca_gwcs = gwcs.WCS([(detector_frame, det2sky), (sky_frame, None)])

        # this has all the fields that fullfovexport is expecting
        asdf.AsdfFile(
            {
                "roman": {
                    "data": image,
                    "dq": dq,
                    "var_rnoise": np.full((4088, 4088), 0.1225, dtype=np.float32),
                    "meta": {
                        "exposure": {"start_time": "2028-12-01T04:00:00"},
                        "ephemeris": {"time": 62106.1666666667},
                        "wcs": sca_gwcs,
                        "instrument": {"optical_element": "F184"},
                    },
                },
                "processinfo": {
                    "medgain": 1.5,
                    "medsky": 0.3,
                    "meta": {
                        "frame_time": frame_time,
                        "read_pattern": read_pattern,
                        "ngrp": ngrp,
                        "K": np.array([0.0, -0.0352734, 0.0, 0.0, 0.0, 0.0352734], dtype=np.float32),
                        "tbar": tbar,
                        "tau": tau,
                    },
                },
            }
        ).write_to(tmp_path + f"/testimage{sca:02d}.asdf")

        if sca not in [1, 8]:
            # Make a mask file except for SCA 1 & 8
            mask = np.zeros((4088, 4088), dtype=np.int16)
            if sca > 10:
                mask[3000:3050, :2044] = -1000
            fits.PrimaryHDU(mask).writeto(tmp_path + f"/testmask{sca:02d}.fits", overwrite=True)

    # make the FFoV image
    FullFoVImage(
        tmp_path + r"/testimage{:02d}.asdf",
        maskfile=tmp_path + r"/testmask{:02d}.fits",
    ).to_file(tmp_path + "/testffov.fits")

    # Now some tests
    with fits.open(tmp_path + "/testffov.fits") as f:
        # tests to run on all the good SCAs
        for sca in range(1, 19):
            if sca != 11:
                d = f[sca].data
                assert 230 <= np.percentile(d, 25) <= 236
                assert 253 <= np.percentile(d, 50) <= 259
                assert 275 <= np.percentile(d, 75) <= 281
                assert 346 <= np.median(d[200:205, 400:800]) <= 353
                for j in range(5):
                    x = np.median(d[100:120, 500 - 25 * j : 520 - 25 * j])
                    if (sca >> j) & 1:
                        assert 1170 <= x <= 1200
                    else:
                        assert x == 1

                # check the masked region
                x = np.median(d[3000:3050, :2044])
                if sca <= 10:
                    assert 253 <= x <= 259
                else:
                    assert x == 0

            # check headers
            h = f[sca].header
            assert h["NAXIS"] == 2
            assert h["NAXIS1"] == 4088
            assert h["NAXIS2"] == 4088
            assert h["EXTNAME"] == f"WFI{sca:02d}"
            if h["ISVALID"]:
                assert h["FILTER"] == "F184"
                assert np.abs(h["MJD"] - 62106.1666666667) < 1.0e-6
                assert 0.45 < h["EQVGAIN"] < 0.46
            else:
                assert sca in [1, 8, 11]  # must be one of the ones we "broke"

        # choose SCA 16 for WCS tests
        sca16wcs = wcs.WCS(f[16].header)
        # lower left
        ra, dec = sca16wcs.pixel_to_world_values(0, 0)
        print(ra, dec)
        assert np.hypot(ra - 149.7268, dec + 19.6893) < 0.001

        # lower right
        ra, dec = sca16wcs.pixel_to_world_values(4087, 0)
        print(ra, dec)
        assert np.hypot(ra - 149.7934, dec + 19.7963) < 0.001

        # upper left
        ra, dec = sca16wcs.pixel_to_world_values(0, 4087)
        print(ra, dec)
        assert np.hypot(ra - 149.6121, dec + 19.7502) < 0.001

        # upper right
        ra, dec = sca16wcs.pixel_to_world_values(4087, 4087)
        print(ra, dec)
        assert np.hypot(ra - 149.6787, dec + 19.8572) < 0.001

    # Similar but check that FullFoVImageFromFile works
    ff = FullFoVImageFromFile(tmp_path + "/testffov.fits")
    valid = np.zeros((18,), dtype=bool)
    for sca in range(1, 19):
        if sca != 11:
            d = ff.hdulist[sca].data
            for j in range(5):
                x = np.median(d[100:120, 500 - 25 * j : 520 - 25 * j])
                if (sca >> j) & 1:
                    assert 1170 <= x <= 1200
                else:
                    assert x == 1
        valid[sca - 1] = ff.hdulist[sca].header["ISVALID"]
    os.remove(tmp_path + "/testffov.fits")

    # Now try again
    valid_target = np.ones((18,), dtype=bool)
    valid_target[11 - 1] = False
    valid_target[8 - 1] = False
    valid_target[1 - 1] = False
    print(valid)
    assert np.all(valid == valid_target)

    # now remove the gain information and check this runs
    with asdf.open(tmp_path + "/testimage05.asdf", mode="rw") as a:
        del a.tree["processinfo"]["medgain"]
        a.update()
    with asdf.open(tmp_path + "/testimage05.asdf") as a:
        assert "medgain" not in a["processinfo"]
    os.remove(tmp_path + "/testffov.fits")
    with pytest.warns(UserWarning, match="Couldn't find median gain, switching to default value."):
        FullFoVImage(
            tmp_path + r"/testimage{:02d}.asdf",
            maskfile=tmp_path + r"/testmask{:02d}.fits",
        ).to_file(tmp_path + "/testffov.fits")

    # Remove old files
    for sca in range(1, 19):
        if sca != 11:
            os.remove(tmp_path + f"/testimage{sca:02d}.asdf")
        if sca not in [1, 8, 11]:
            os.remove(tmp_path + f"/testmask{sca:02d}.fits")
