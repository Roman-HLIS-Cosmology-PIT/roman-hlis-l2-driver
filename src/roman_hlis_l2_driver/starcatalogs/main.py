"""Driver for masking bright star spikes."""


import asdf
import numpy as np
from pyimcom.wcsutil import PyIMCOM_WCS

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

    # set up coordinates
    _s = np.array(range(nside)).astype(np.float32)
    dx, dy = np.meshgrid(_s, _s)
    dx -= x
    dy -= y
    slope = np.float32(np.tan(wedge))

    for j in range(n):
        rotx = dx * np.cos(sp_theta[j]) + dy * np.sin(sp_theta[j])
        roty = np.abs(dy * np.cos(sp_theta[j]) - dx * np.sin(sp_theta[j]))
        mask |= (
            (rotx > 0)
            & (rotx < sp_length[j] + 0.5)
            & (roty < sp_width[j] + 0.5)
            & (roty < slope * rotx - 0.5)
        )

    return mask


def stardata_to_mask(use_filter, sca, x, y, thresh):
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

    Returns
    -------
    np.ndarray of bool
        An array of shape (`nside`, `nside`); values are True if masked.

    """

    n = len(x)
    nside = 4088

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
