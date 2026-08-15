"""
Tools to build a star catalog internally from the images.

"""

import glob
import re

import asdf
import numpy as np
import numpy.lib.recfunctions as rfn
import sep
from roman_datamodels.dqflags import pixel


def encirc_center(arr):
    """
    Finds the center and radius of the smallest circle enclosing the "True" values in an array.

    Parameters
    ----------
    arr : np.ndarray of bool or bool-like
        The input array. Must be 2 dimensional.

    Returns
    -------
    x, y : float
        The position of the circumcenter (zero offset).
    r : float
        The circumradius.

    Notes
    -----
    If the array is all False, returns -1., -1., -1.

    """

    # Get the identified points
    ny, nx = np.shape(arr)
    count = np.count_nonzero(arr)
    yf, xf = np.where(arr)
    x_, y_ = np.meshgrid(np.arange(nx), np.arange(ny))

    # first handle contingencies for <3 points
    if count == 0:
        return -1.0, -1.0, -1.0
    if count == 1:
        return float(xf[0]), float(yf[0]), 0.0
    if count == 2:
        r = np.hypot(xf[0] - xf[1], yf[0] - yf[1]) / 2.0
        return 0.5 * (xf[0] + xf[1]), 0.5 * (yf[0] + yf[1]), r

    # Now find the circumscribing points in order.
    # First get "top right" in ds9.
    pos = np.argmax(np.arange(ny * nx) * arr.ravel())
    yp = [pos // nx]
    xp = [pos % nx]
    nc = (ny + nx) ** 2  # max for rank-ordering
    searchdir = (0, -1)  # ordering: dy, dx. *must be integer*
    closed = False
    while not closed:
        gr = searchdir[0] * (y_ - yp[-1]) + searchdir[1] * (x_ - xp[-1])
        gi = searchdir[1] * (y_ - yp[-1]) - searchdir[0] * (x_ - xp[-1])
        rank = -np.angle(-(gr + 1j * gi) * np.exp(1j / nc)) + np.hypot(x_ - xp[-1], y_ - yp[-1]) / nc**2
        rank[yp[-1], xp[-1]] = -1000
        pos2 = np.argmax(np.where(arr.ravel(), rank.ravel(), -1000))
        yp += [pos2 // nx]
        xp += [pos2 % nx]
        dy = yp[-1] - yp[-2]
        dx = xp[-1] - xp[-2]
        g = np.gcd(dy, dx)
        searchdir = (dy // g, dx // g)
        if xp[0] == xp[-1] and yp[0] == yp[-1]:
            closed = True
        if len(xp) == 16384:
            raise ValueError("Polygon convergence failed.")
    npts = len(xp) - 1  # number of points in the cicrumscribing polygon

    # now get the candidate centers
    nbis = npts * (npts - 1) // 2
    ntrip = npts * (npts - 1) * (npts - 2) // 6
    ntot = nbis + ntrip
    pos = np.zeros((ntot, 2))
    ipos = 0
    for i in range(1, npts):
        for j in range(i):
            pos[ipos, :] = s = np.array([(yp[i] + yp[j]) / 2.0, (xp[i] + xp[j]) / 2.0])
            p = np.array([-(xp[i] - xp[j]) / 2.0, (yp[i] - yp[j]) / 2.0])
            ipos += 1
            for k in range(j):
                # now we go from s a distance d along p until we are the same
                # distance from point i and point k.
                #
                # formula: dist[from i]^2 - dist[from k]^2 = |s + d*p - r_i|^2 - |s + d*p - r_k|^2
                # = 2 (s + d*p - (r_i+r_k)/2.) dot (r_k - r_i).
                #
                # setting this equal to zero gives
                # d = [(r_k - r_j) dot (r_k - r_i)] / [2 p dot (r_k - r_i)]
                rki = np.array([yp[k] - yp[i], xp[k] - xp[i]])
                rkj = np.array([yp[k] - yp[j], xp[k] - xp[j]])
                d = np.dot(rkj, rki) / 2.0 / np.dot(p, rki)
                pos[ipos, :] = s + d * p
                ipos += 1

    # now test the radii
    rtest = np.zeros((ntot,))
    for ipos in range(ntot):
        rtest[ipos] = np.amax(np.hypot(y_ - pos[ipos, 0], x_ - pos[ipos, 1]) * arr)
    ipos = np.argmin(rtest)
    return pos[ipos, 1], pos[ipos, 0], rtest[ipos]


def brightobj_from_scaimg(infile, obsid=-1, thresh=50.0, verbose=False):
    """
    Extract bright objects from an SCA image.

    Parameters
    ----------
    infile : str or str-like
        The input image.
    obsid : int, optional
        The observation ID (defaults to -1).
    thresh : float, optional
        The threshold to use (in DN/s).
    verbose : bool, optional
        Whether to print the object catalogs and data to the terminal.

    Returns
    -------
    np.recarray
        Bright object catalog from ``sep.extract``. The most important fields are:
        * ``ra``, ``dec``: The object position (in degrees).
        * ``xcc``, ``ycc``: The center as determined from the isophote.
        * ``eflux``: The estimated flux (in total DN/s).

    Notes
    -----
    Don't use these estimates for good astrometry or photometry! They are intended
    only for masking or initial guesses, and errors of a few tenths of a pixel or
    tens of percents in photometry are common. Remember, this routine doesn't fit the object
    itself, instead it is using the wings so that it can get an estimate even for objects
    that are completely saturated.

    """

    maskbits = (
        pixel.JUMP_DET
        | pixel.DROPOUT
        | pixel.GW_AFFECTED_DATA
        | pixel.AD_FLOOR
        | pixel.NON_SCIENCE
        | pixel.DEAD
        | pixel.HOT
        | pixel.LOW_QE
        | pixel.TELEGRAPH
        | pixel.NO_FLAT_FIELD
        | pixel.NO_GAIN_VALUE
        | pixel.NO_SAT_CHECK
        | pixel.UNRELIABLE_BIAS
        | pixel.UNRELIABLE_DARK
        | pixel.UNRELIABLE_SLOPE
        | pixel.UNRELIABLE_FLAT
        | pixel.UNRELIABLE_RESET
        | pixel.OTHER_BAD_PIXEL
    )
    with asdf.open(infile) as f:
        img = np.array(f["roman"]["data"])
        (ny, nx) = np.shape(img)
        es = np.array(f["processinfo"]["endslice"])
        mask = np.array(f["roman"]["dq"]) & maskbits != 0
        # pixels flagged for saturating immediately are OK for this
        mask |= np.logical_and(np.array(f["roman"]["dq"]) & pixel.DO_NOT_USE != 0, es != 0)
        element = f["roman"]["meta"]["instrument"]["optical_element"]
        wcs = f["roman"]["meta"]["wcs"]
        sca = int(f["roman"]["meta"]["instrument"]["detector"][-2:])

    wl = float(element[-3:]) / 100.0  # wavelength in microns
    ldp = wl * 0.79  # lambda/D in pixels
    alpha = ldp * 2.0 / np.pi / (1 - 0.31)  # 0.31 = obscuration
    # this is an effective Moffat width for beta=3/2
    # that gets the right diffraction ring settings
    if verbose:
        print("alpha =", alpha, "pixels")

    # re-center and grow object if necessary
    obj, map = sep.extract(img, thresh, mask=mask, deblend_cont=0.5, segmentation_map=True)
    N = np.shape(obj)[0]
    xcc = np.zeros(N)
    ycc = np.zeros(N)
    rcc = np.zeros(N)
    eflux = np.zeros(N)
    idsca = np.full(N, 100 * obsid + sca, dtype=np.int64)
    xi, yi = np.meshgrid(np.arange(nx), np.arange(ny))
    for i in range(N):
        oi = obj[i]
        x_, y_, rcc[i] = encirc_center(map[oi["ymin"] : oi["ymax"] + 1, oi["xmin"] : oi["xmax"] + 1] == i + 1)
        xcc[i] = x_ + oi["xmin"]
        ycc[i] = y_ + oi["ymin"]
        r = np.hypot(xi - xcc[i], yi - ycc[i])
        fluxfrac = 1.0 / (2 * np.pi * alpha) * (1.0 + r**2 / alpha**2) ** (-1.5)
        ring = np.where(
            np.logical_and(
                np.logical_and(r > np.pi / 4.0 * rcc[i], r < np.pi / 2.0 * rcc[i]), np.logical_not(mask)
            )
        )
        if len(ring[0]) == 0:
            ring = np.where(
                np.logical_and(
                    np.logical_and(r > np.pi / 4.0 * rcc[i], r < np.pi / 2.0 * rcc[i]), map == i + 1
                )
            )
        if len(ring[0]) == 0:
            # This isn't supposed to happen.
            print("ERROR")
            print(map[oi["ymin"] : oi["ymax"] + 1, oi["xmin"] : oi["xmax"] + 1] == i + 1)
            print("x range", oi["xmin"], oi["xmax"], "y range", oi["ymin"], oi["ymax"])
            print((xcc[i], ycc[i]), "r =", rcc[i])
            print("Mask:")
            print(mask[oi["ymin"] : oi["ymax"] + 1, oi["xmin"] : oi["xmax"] + 1])
            print("---")
            raise ValueError("Invalid mask.")
        eflux[i] = max(np.mean(img[ring]) / np.mean(fluxfrac[ring]), obj[i]["flux"])
    obj = rfn.append_fields(
        obj,
        ["idsca", "xcc", "ycc", "rcc", "eflux"],
        [idsca, xcc, ycc, rcc, eflux],
        usemask=False,
        asrecarray=True,
    )

    # now get the RA & Dec
    ra, dec = wcs.pixel_to_world_values(obj["xcc"], obj["ycc"])
    obj = rfn.append_fields(obj, ["ra", "dec"], [ra, dec], usemask=False, asrecarray=True)

    if verbose:
        for i in range(N):
            print(
                f"{obj[i]["xcc"]:7.2f} {obj[i]["ycc"]:7.2f} {obj[i]["rcc"]:7.2f} "
                f"{obj[i]["flux"]:9.1f} {obj[i]["eflux"]:9.1f}"
            )

    # this was for testing only
    # maskimg = fits.ImageHDU(np.where(~mask, img, -100.0))
    # fits.HDUList([fits.PrimaryHDU(map), fits.ImageHDU(es), maskimg]).writeto("test.fits", overwrite=True)

    return obj


def cluster(objcat, dist, verbose=False):
    """
    Constructs clusters from a table of objects.

    The table is assumed to be pre-sorted from lowest to highest priority.

    Parameters
    ----------
    objcat : np.recarray
        A record array, including "ra" and "dec" fields (in degrees).
    dist : float
        Linking distance (in degrees).
    verbose : bool, optional
        Print the matches as they are found.

    Returns
    -------
    np.ndarray of int
        A vector of integers of the same length as objcat, indicating
        the parent that the object is associated with. A "-1" indicates
        that there is no parent.
    np.ndarray of int
        A vector of integers of the same length as objcat, indicating
        the number of child objects.

    Notes
    -----
    This function is not super efficient; it uses an N^2 algorithm, since
    its intended use in the pipeline is on a small number of sources, and
    its total cost is negligible.

    """

    # Get 3D coordinates and 3D distance
    deg = np.pi / 180.0
    N = len(objcat)
    pos = np.zeros((N, 3))
    pos[:, 0] = np.cos(objcat["dec"] * deg) * np.cos(objcat["ra"] * deg)
    pos[:, 1] = np.cos(objcat["dec"] * deg) * np.sin(objcat["ra"] * deg)
    pos[:, 2] = np.sin(objcat["dec"] * deg)
    r = 2 * np.sin(dist * deg / 2.0)

    # build output array
    parent = np.full(N, -1, dtype=np.int32)
    nchild = np.zeros(N, dtype=np.int32)

    # Search for children of object k, in order from brightest to faintest.
    # We skip the faintest, since it can't have any children (hence
    # the for loop has N-1, 0, -1 instead of N-1, -1, -1).
    for k in range(N - 1, 0, -1):
        # the parent to assign to
        assign_parent = k if parent[k] == -1 else parent[k]

        # get the matches
        sep = np.linalg.norm(pos[:k, :] - pos[k, :][None, :], axis=1)
        match = np.where(sep < r)[0]
        if len(match) > 0:
            if verbose:
                print(k, match, np.arcsin(sep[match] / 2) * 2 / deg)
            nchild[assign_parent] += len(match)
            parent[match] = assign_parent

    return parent, nchild


def brightobj_from_manyimg(
    infile_format, thresh=50.0, maximg=None, matchrad=0.0002777777777777778, verbose=False, clean=False
):
    """
    Extract bright objects from a group of images.

    Parameters
    ----------
    infile_format : str
        The input file as a formatted string. One should be able to write
        ``infformat.format(obsid, sca)``
    thresh : float, optional
        The threshold to use (in DN/s).
    maximg : int, optional
        Use a maximum of this many SCAs; primarily used for testing.
    matchrad : float, optional
        Matching radius for the catalog (in degrees).
    verbose : bool, optional
        Whether to print the object catalogs and data to the terminal.
    clean : bool, optional
        Whether to remove unmatched objects (e.g., only one observation)
        from the catalog.

    Returns
    -------
    np.recarray
        Bright object catalog from ``sep.extract``. The most important fields are:
        * ``idsca``: The object ID/SCA of the detection (in the form ``100*obsid+sca``).
        * ``ra``, ``dec``: The object position (in degrees).
        * ``eflux``: The estimated flux (in total DN/s).

    See Also
    --------
    brightobj_from_scaimg
        Similar, but for only a single image.

    """

    # take apart the format
    matchstr = re.sub(r"{[^{}]+}", r"({[^{}]+})", infile_format)
    ingrps = re.match(matchstr, infile_format).groups()
    ind = []
    ngrp = len(ingrps)
    if ngrp < 2:
        raise ValueError("Can't find obsid and sca in string format: " + infile_format)
    for j in range(ngrp):
        s = re.match(r"{(\d*):", ingrps[j])
        if not s:
            raise ValueError("Can't find obsid and sca in string format: " + infile_format)
        s = s.group(1)
        if s == "":
            ind.append(j)
        else:
            ind.append(int(s))
    ind_obsid = -1
    ind_sca = -1
    for j in range(ngrp):
        if ind[j] == 0:
            ind_obsid = j
        if ind[j] == 1:
            ind_sca = j
    if ind_obsid < 0 or ind_sca < 0:
        raise ValueError("Can't find obsid and sca in string format: " + infile_format)

    # now get the files
    matchstr = re.sub(r"{[^{}]+}", r"(\\d+)", infile_format)
    filelist_coarse = glob.glob(re.sub(r"{[^{}]+}", r"*", infile_format))
    idsca = []
    for fn in filelist_coarse:
        m = re.match(matchstr, fn)
        if m:
            idsca.append((int(m.groups()[ind_obsid]), int(m.groups()[ind_sca])))
    idsca = idsca[:maximg]  # truncate if needed
    if verbose:
        print("Searching:", idsca)

    _is_init = False
    for j in range(len(idsca)):
        fn = infile_format.format(*idsca[j])
        cat = brightobj_from_scaimg(fn, obsid=idsca[j][0], thresh=thresh, verbose=verbose)
        if len(cat) > 1:
            if _is_init:
                ccat = np.concatenate((ccat, cat))  # ccat will be defined before we get here # noqa: F821
            else:
                ccat = cat
                _is_init = True
        if verbose:
            print(fn)
    ccat.sort(order="eflux")

    parent, nchild = cluster(ccat, matchrad, verbose=verbose)
    ccat = rfn.append_fields(
        ccat,
        ["objid", "parent", "nchild"],
        [np.arange(len(ccat)), parent, nchild],
        usemask=False,
        asrecarray=True,
    )

    if clean:
        ccat = ccat[np.where(np.logical_or(ccat["parent"] != -1, ccat["nchild"] > 0))]
    return ccat
