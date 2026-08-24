"""Driver for masking bright star spikes."""


import sys

import asdf
import numpy as np
from pyimcom.wcsutil import PyIMCOM_WCS

from .internal import brightobj_from_manyimg
from .interp import extract_spikedata


def spikedata_to_mask(nside, x, y, sp_theta, sp_length, sp_width, wedge):
    """
    Makes a spike mask.

    Parameters
    ----------
    nside : int
        The side length of the SCA.
    x, y : float
        The location of the bright star to mask.
    sp_theta, sp_length, sp_width : np.ndarray of float
        Spike angles, lengths, and half-widths; shape (n,).
    wedge : float
        The half-opening angle to consider for the spikes (so we don't mask
        the center where they all overlap).

    Returns
    -------
    np.ndarray of bool
        An array of shape (`nside`, `nside`); values are True if masked.

    """

    mask = np.zeros((nside, nside), dtype=bool)
    n = len(sp_theta)

    # replace with short versions
    sp_theta = sp_theta.astype(np.float32)
    sp_length = sp_length.astype(np.float32)
    sp_width = sp_width.astype(np.float32)

    # get bounding boxes
    hyp = np.amax(np.hypot(sp_length, sp_width))
    xmin = max(int(np.floor(x - hyp)), 0)
    xmax = min(int(np.floor(x + hyp + 1)), nside)
    ymin = max(int(np.floor(y - hyp)), 0)
    ymax = min(int(np.floor(y + hyp + 1)), nside)

    # only execute if bounding box not empty
    if xmax > xmin and ymax > ymin:
        # set up coordinates
        _sx = np.array(range(xmin, xmax)).astype(np.float32)
        _sy = np.array(range(ymin, ymax)).astype(np.float32)
        dx, dy = np.meshgrid(_sx, _sy)
        dx -= x
        dy -= y
        slope = np.float32(np.tan(wedge))

        for j in range(n):
            rotx = dx * np.cos(sp_theta[j]) + dy * np.sin(sp_theta[j])
            roty = np.abs(dy * np.cos(sp_theta[j]) - dx * np.sin(sp_theta[j]))
            mask[ymin:ymax, xmin:xmax] |= (
                (rotx > 0)
                & (rotx < sp_length[j] + 0.5)
                & (roty < sp_width[j] + 0.5)
                & (roty < slope * rotx - 0.5)
            )

    return mask


def stardata_to_mask(use_filter, sca, x, y, thresh, nside=4088):
    """
    Makes a spike mask for a star catalog in Science coordinates.

    Parameters
    ----------
    use_filter : char
        The filter character code.
    sca : int
        The SCA number (1..18).
    x, y : np.ndarray of float
        The location of the bright star to mask.
    thresh : np.ndarray of float
        The threshold (relative to star brightness) to mask.
    nside : int, optional
        Side length of the SCA (4088 for Roman).

    Returns
    -------
    np.ndarray of bool
        An array of shape (`nside`, `nside`); values are True if masked.

    """

    n = len(x)

    mask = np.zeros((nside, nside), dtype=bool)
    for j in range(n):
        mask |= spikedata_to_mask(
            nside, x[j], y[j], *extract_spikedata(use_filter, sca, x[j], y[j], thresh[j])
        )

    return mask


def stardata_to_mask_wcs(l2file, catalog, thresh, dp=256.0):
    """
    Makes a spike mask for a star catalog in world coordinates.

    Parameters
    ----------
    l2file : str
        The input ASDF L2 file for this exposure.
    catalog : np.recarray
        The star catalog; should have at least ra, dec, and estar.
    thresh : np.ndarray of float
        The threshold of diffraction spike brightness to mask.
    dp : float, optional
        The distance in pixels around the image to consider for masking stars.

    Returns
    -------
    np.ndarray of bool
        An array of shape (`nside`, `nside`); values are True if masked.

    """

    filterkeys = {
        "F062": "R",
        "F087": "Z",
        "F106": "Y",
        "F129": "J",
        "F158": "H",
        "F184": "F",
        "F213": "K",
        "F146": "W",
    }

    # extract information from the file
    with asdf.open(l2file) as a:
        use_filter = filterkeys[a["roman"]["meta"]["instrument"]["optical_element"]]
        sca = int(a["roman"]["meta"]["instrument"]["detector"][-2:])
        pwcs = PyIMCOM_WCS(a["roman"]["meta"]["wcs"])
        nside = np.shape(a["roman"]["data"])[-1]

    assert nside == 4088  # REMOVE

    # extract only stars near that SCA
    ctr = pwcs.all_pix2world([[2043.5, 2043.5]], 0)[0]
    ra_ctr = ctr[0]
    dec_ctr = ctr[1]
    x_ctr = np.cos(np.radians(dec_ctr)) * np.cos(np.radians(ra_ctr))
    y_ctr = np.cos(np.radians(dec_ctr)) * np.sin(np.radians(ra_ctr))
    z_ctr = np.sin(np.radians(dec_ctr))
    x = np.cos(np.radians(catalog["dec"])) * np.cos(np.radians(catalog["ra"]))
    y = np.cos(np.radians(catalog["dec"])) * np.sin(np.radians(catalog["ra"]))
    z = np.sin(np.radians(catalog["dec"]))
    sep = np.sqrt((x - x_ctr) ** 2 + (y - y_ctr) ** 2 + (z - z_ctr) ** 2)
    sep = np.arcsin(sep / 2.0) * 2.0
    subcatalog = catalog[sep < np.radians(0.1)]  # 0.1 degree radius

    x, y = pwcs.all_world2pix(subcatalog["ra"], subcatalog["dec"], 0)
    trim = np.logical_and(np.abs(x - 2043.5) < 2044.5 + dp, np.abs(y - 2043.5) < 2044.5 + dp)

    return stardata_to_mask(use_filter, sca, x[trim], y[trim], thresh / subcatalog["eflux"][trim])


def stardata_to_mask_manyfiles(l2files, catalog, thresh, dp=256.0, update=True):
    """
    Makes a spike mask for a star catalog in world coordinates.

    Parameters
    ----------
    l2files : list of str
        The input ASDF L2 file for this exposure.
    catalog : np.recarray
        The star catalog; should have at least ra, dec, and estar.
    thresh : float
        The threshold of diffraction spike brightness to mask.
    dp : float, optional
        The distance in pixels around the image to consider for masking stars.
    update : bool, optional
        Whether to update the "mask" 0x8 bit in the ASDF files. (You should only need to
        turn this off for diagnostics.)

    Returns
    -------
    np.ndarray of bool
        Table of the masked pixels in each exposure; shape (nexp, nside, nside).

    """

    # Extract setup information
    nexp = len(l2files)
    with asdf.open(l2files[0]) as a:
        nside = np.shape(a["roman"]["data"])[-1]
    mask = np.zeros((nexp, nside, nside), dtype=bool)

    for j in range(nexp):
        mask[j, :, :] = stardata_to_mask_wcs(l2files[j], catalog, thresh, dp=dp)
        print(f"{j}/{nexp}: {l2files[j]} {np.count_nonzero(mask[j,:,:])}")
        sys.stdout.flush()

    if update:
        for j in range(nexp):
            with asdf.open(l2files[j], mode="rw") as a:
                a["mask"] &= ~np.uint8(0x8)
                a["mask"] |= mask[j].astype(np.uint8) << 3
                a["processinfo"]["spike_complete"] = True
                a["processinfo"]["spike_pars"] = {"thresh": thresh, "dp": dp}
                a.update()

    return mask


def spike_driver(
    infile_format,
    thresh_mask,
    dp=256.0,
    update=True,
    thresh_detect=50.0,
    maximg=None,
    matchrad=0.0002777777777777778,
):
    """
    Extract bright objects from a group of images.

    Parameters
    ----------
    infile_format : str
        The input file as a formatted string. One should be able to write
        ``infformat.format(obsid, sca)``
    thresh_mask : float
        The threshold of diffraction spike brightness to mask.
    dp : float, optional
        The distance in pixels around the image to consider for masking stars.
    update : bool, optional
        Whether to update the "mask" 0x8 bit in the ASDF files. (You should only need to
        turn this off for diagnostics.)
    thresh_detect : float, optional
        The threshold to use (in DN/s).
    maximg : int, optional
        Use a maximum of this many SCAs; primarily used for testing.
    matchrad : float, optional
        Matching radius for the catalog (in degrees).

    Returns
    -------
    stars_unique : np.recarray
        Bright object catalog from ``sep.extract``. The most important fields are:
        * ``idsca``: The object ID/SCA of the detection (in the form ``100*obsid+sca``).
        * ``ra``, ``dec``: The object position (in degrees).
        * ``eflux``: The estimated flux (in total DN/s).
    npix_mask : int
        The total number of pixels masked.

    """

    stars_recovered, idsca = brightobj_from_manyimg(infile_format, thresh=thresh_detect, clean=True)
    stars_unique = stars_recovered[stars_recovered["nchild"] > 0]

    l2files = [infile_format.format(*i) for i in idsca]
    print(l2files)
    m = stardata_to_mask_manyfiles(l2files, stars_unique, thresh_mask, dp=dp, update=update)
    npix_mask = np.count_nonzero(m)

    return stars_unique, npix_mask
