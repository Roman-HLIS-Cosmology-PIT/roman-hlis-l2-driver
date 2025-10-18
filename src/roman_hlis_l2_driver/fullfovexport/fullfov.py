"""
Full field of view handling tools

Classes
-------
FullFoVImage
    An image of the full field of view.
FullFoVImageFromFile
    Same as FullFoVImage, but loaded from disk.

"""

import asdf
import numpy as np
from astropy.io import fits
from pyimcom.wcsutil import LocWCS

from .. import pars


def get_t_eff(K, tau, tbar):
    """
    Compute effective exposure time for multiple reads.

    Parameters
    ----------
    K : array-like
        The weights of each resultant as a 1D array.
    tau : array-like
        The variance-weighted time-stamp of each resultant as a 1D array.
    tbar : array-like
        The mean time-stamp of each resultant as a 1D array.

    Returns
    -------
    float
        The effective time, given by <I>/Var(I) where I is the slope in e/s.

    Notes
    -----
    The `tau` and `tbar` variables are as defined in Casertano,
    "Determining the best-fitting slope and its uncertainty
    for up-the-ramp sampled images
    with unevenly distributed resultants" (2022),
    Technical Report Roman-STScI-000394.

    The `K` array has units of 1/s, and is used in the form of
    I = K dot Q, where Q is the vector of 1D collected charge at
    each time stamp.

    """

    # re-format as arrays
    K = np.array(K)
    tau = np.array(tau)
    tbar = np.array(tbar)

    # get covariance matrix for current I=1.
    ngrp = len(K)
    C = np.zeros(ngrp, ngrp)
    for i in range(ngrp):
        C[i, :i] = C[:i, i] = tbar[:i]
        C[i, i] = tau[i]
    # and the variance for the compressed slope
    var = np.dot(K, C @ K)

    return 1.0 / var


class FullFoVImage:
    """
    An image of the full field of view.

    Parameters
    ----------
    infile : str
        A formatted string, so that ``infile.format(sca)`` returns
        the path of the SCA (with sca as an integer starting from 1).
    maskfile : str, optional
        A formatted string, so that ``maskfile.format(sca)`` returns
        the path of the SCA (with sca as an integer starting from 1).
        If not provided, then no mask is applied (not recommended).
    softbias : int, optional
        The bias (code corresponding to zero count rate).
    slopemax : float, optional
        The slope in linearized DN/s corresponding to 16-bit saturation.
        (Default is likely good for the HLWAS.)
    wcs_src : str, optional
        Source to take WCS (default is from `infile`).

    Attributes
    ----------
    hdulist : astropy.io.fits.HDUList
        The FITS HDUList representation of the images + metadata.

    Methods
    -------
    to_file
        Saves the images to a FITS file.

    """

    def __init__(self, infile, maskfile=None, softbias=256, slopemax=700.0, wcs_src="L2"):
        # General information goes in the Primary HDU header.
        phdu = fits.PrimaryHDU()
        phdu.header["SOFTBIAS"] = (softbias, "Code corresponding to sky level.")
        dslope = slopemax / (2**16 - 2 - softbias)
        slopemin = slopemax - dslope * (2**16 - 3)
        phdu.header["SLOPEMIN"] = (slopemin, "Slope in DN_lin/s corresponding to 1.")
        phdu.header["SLOPEMAX"] = (slopemax, "Slope in DN_lin/s corresponding to 65534.")
        phdu.header["DSLOPE"] = (dslope, "Slope unit in DN_lin/s.")
        hdulist = [phdu]

        n_side = pars.n_side  # this is just for convenience since we'll use it so much

        for sca in range(1, 1 + pars.n_sca):
            im = np.zeros((n_side, n_side), dtype=np.uint16)
            valid = True
            self.wcs = None
            try:
                with asdf.open(infile.format(sca)) as a:
                    im[:, :] = np.clip(np.rint(softbias + a["roman"]["data"] / dslope), 1, 2**16 - 1).astype(
                        np.uint16
                    )
                    start_time = a["roman"]["meta"]["exposure"]["start_time"]
                    meta = a["roman"]["meta"]  # for shorthand
                    mjd = meta["ephemeris"]["time"]

                    # effective gain information
                    medgain = a["processinfo"]["medgain"]
                    t_eff = get_t_eff(meta["K"], meta["tbar"], meta["tau"])
                    print("t_eff =", t_eff, "s")

                    # the WCS
                    if wcs_src.lower() == "l2":
                        self.wcs = LocWCS(meta["wcs"], N=n_side)
                        self.wcs.wcs_approx_sip(p_order=4)

            except FileNotFoundError:
                valid = False
            mask = False
            if maskfile is not None:
                try:
                    with fits.open(maskfile.format(sca)) as m:
                        im[:, :] = np.where(m[0].data >= -999, im, 0)
                    mask = True
                except FileNotFoundError:
                    valid = False
                # masked data is now 0's

            new_hdu = fits.ImageHDU(im, name=f"WFI{sca:02d}")
            new_hdu.header["ISVALID"] = (valid, "Was this SCA found?")
            new_hdu.header["HASMASK"] = (mask, "Was a mask applied?")
            new_hdu.header["HASWCS"] = bool(self.wcs is not None)

            # for now, not using the pixel-level error map
            new_hdu.header["ERRMAP"] = ("NULL", "Error map name")

            # add noise information
            eqvgain = dslope * t_eff * medgain
            bkgndvar1 = np.median(a["roman"]["var_rnoise"]) / dslope**2
            bkgndvar2 = a["processinfo"]["medsky"] / (t_eff * medgain * dslope**2)
            print(f"g_eqv = {eqvgain}, variance {bkgndvar1} from read, {bkgndvar2} from sky Poisson")
            new_hdu.header["EQVGAIN"] = (eqvgain, "Equivalent gain in weighted electrons per int in file")
            new_hdu.header["BKGNDVAR"] = (bkgndvar1 + bkgndvar2, "Variance at background level")

            # add WCS
            if self.wcs is not None:
                new_hdu.header.update(self.wcs.approx_wcs.to_header(relax=True))
                new_hdu.header["MAXWCSER"] = (self.wcs.wcs_max_err, "max of error map in pixels")
            hdulist.append(new_hdu)

        # save collected metadata
        phdu.header["MJD"] = mjd
        phdu.header["TSTART"] = str(start_time)

        self.hdulist = fits.HDUList(hdulist)  # save this as an HDUList

    def to_file(self, filename, overwrite=False):
        """
        Saves to a FITS file.

        Parameters
        ----------
        filename : str or str-like
            The file name to save the FITS images.
        overwrite : bool, optional
            Whether to overwrite the destination file.

        """

        self.hdulist.writeto(filename, overwrite=overwrite)


class FullFoVImageFromFile(FullFoVImage):
    """
    Like `FullFovImage`, but construct from a saved file instead of re-doing the computations.

    Parameters
    ----------
    infile : str or str-like
        A valid FITS file containing the full field of view.

    See Also
    --------
    FullFoVImage
        Base class.

    """

    def __init__(self, infile):
        with fits.open(infile, memmap=False) as f:
            self.hdulist = f
            # turn off memmap so that there are no references to disk after this
