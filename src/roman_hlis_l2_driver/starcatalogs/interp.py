"""Interpolation routines specific to the stars code."""


from importlib.resources import files

import asdf
import numpy as np


def interp5pt(arr, u, v, axis):
    """
    Interpolates an array of length 5 to a point (u, v).

    Parameters
    ----------
    arr : np.ndarray
        The input array.
    u, v : float
        The location to interpolate to. Inputs are assumed to be at:

        - ``arr[...,0,...]`` = value at (u=-1, v=-1)
        - ``arr[...,1,...]`` = value at (u=1, v=-1)
        - ``arr[...,2,...]`` = value at (u=0, v=0)
        - ``arr[...,3,...]`` = value at (u=-1, v=1)
        - ``arr[...,4,...]`` = value at (u=1, v=1)

    axis : int
        Which axis to compress along.

    Returns
    -------
    np.ndarray
        Has the indicated axis of `arr` squashed.

    """

    # interpolation weights - quadratic from each point
    weights = np.array(
        [
            (-u - v + u * v) / 4.0 + 0.125 * (u**2 + v**2),
            (u - v - u * v) / 4.0 + 0.125 * (u**2 + v**2),
            1.0 - 0.5 * (u**2 + v**2),
            (-u + v - u * v) / 4.0 + 0.125 * (u**2 + v**2),
            (u + v + u * v) / 4.0 + 0.125 * (u**2 + v**2),
        ]
    )

    return np.tensordot(arr, weights, axes=([axis], [0]))


def extract_data(spikedata, use_filter, sca, x, y, thresh):
    """
    Interpolates spike data from a dictionary.

    Parameters
    ----------
    spikedata : dict or dict-like
        Contains the attributes:

        - filters : str
              The string of filter characters in the order of the table.
        - dsp : float
              How far off the edge of the chip we have placed the interpolating points.
        - wedge : float
              The half-opening angle within which we consider each spike.
        - thresh : np.ndarray of float
              The pre-tabulated thresholds, in decreasing order.
        - angles : np.ndarray of float
              The angles of the spikes.
        - lengths, widths : np.ndarray of float
              The angles of the spikes.
    use_filter : char
        The filter to use; one of "RZYJHFKW".
    sca : int
        The SCA number (1..18).
    x, y : float
        Science coordinate positions.
    thresh : float
        The masking threshold to use (interpolates from the table).

    Returns
    -------
    angles, lengths, widths : np.ndarray of float
        The angles, lengths, and widths of the spikes to mask.
        Each an array of length 12.
    wedge : float
        The half-opening angle within which we consider each spike (pass-through).

    """

    # Where to interpolate to?
    u = (x - 2043.5) / (2044.5 + spikedata["dsp"])
    v = (y - 2043.5) / (2044.5 + spikedata["dsp"])

    # weight index
    th = spikedata["thresh"]
    w = np.zeros(len(th))
    if thresh >= th[0]:
        w[0] = 1.0
    elif thresh <= th[-1]:
        # this extrapolates beyond the PSFSim model grid with a 1/3 power law,
        # which is appropriate for diffraction from a slightly curved edge.
        w[-1] = (th[-1] / thresh) ** (1 / 3)
    else:
        j = np.amax(np.where(th > thresh))
        fr = np.log(thresh / th[j]) / np.log(th[j + 1] / th[j])
        w[j] = 1 - fr
        w[j + 1] = fr

    # Now extract the interpolated data
    ifilt = spikedata["filters"].index(use_filter)
    angles = interp5pt(spikedata["angles"][ifilt, sca - 1, ...], u, v, 0)
    lengths = interp5pt(spikedata["lengths"][ifilt, sca - 1, ...], u, v, 0) @ w
    widths = interp5pt(spikedata["widths"][ifilt, sca - 1, ...], u, v, 0) @ w

    return angles, lengths, widths, spikedata["wedge"]


def extract_spikedata(use_filter, sca, x, y, thresh):
    """
    Interpolates spike data from stored data files.

    Parameters
    ----------
    use_filter : char
        The filter to use; one of "RZYJHFKW".
    sca : int
        The SCA number (1..18).
    x, y : float
        Science coordinate positions.
    thresh : float
        The masking threshold to use (interpolates from the table).

    Returns
    -------
    angles, lengths, widths : np.ndarray of float
        The angles, lengths, and widths of the spikes to mask.
        Each an array of length 12.

    """

    with asdf.open(files("roman_hlis_l2_driver.data").joinpath("spikedata.asdf")) as spikedata:
        return extract_data(spikedata, use_filter, sca, x, y, thresh)
