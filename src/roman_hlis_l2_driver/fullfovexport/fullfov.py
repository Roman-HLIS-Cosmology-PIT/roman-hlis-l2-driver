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

from .. import pars


class FullFoVImage:
    """
    An image of the full field of view.

    Parameters
    ----------
    infile : str
        A formatted string, so that ``infile.format(sca)`` returns
        the path of the SCA (with sca as an integer starting from 1).
    softbias : int, optional
        The bias (code corresponding to zero count rate).
    slopemax : float, optional
        The slope in linearized DN/s corresponding to 16-bit saturation.
        (Default is likely good for the HLWAS.)

    Attributes
    ----------
    hdulist : astropy.io.fits.HDUList
        The FITS HDUList representation of the images + metadata.

    Methods
    -------
    to_file
        Saves the images to a FITS file.

    """

    def __init__(self, infile, softbias=256, slopemax=700.0):
        # General information goes in the Primary HDU header.
        phdu = fits.PrimaryHDU()
        phdu.header["SOFTBIAS"] = (softbias, "Signal corresponding to sky level.")
        dslope = slopemax / (2**16 - 1 - softbias)
        slopemin = slopemax - dslope * (2**16 - 1)
        phdu.header["SLOPEMIN"] = (slopemin, "Slope in DN_lin/s corresponding to 0.")
        phdu.header["SLOPEMAX"] = (slopemax, "Slope in DN_lin/s corresponding to 65535.")
        hdulist = [phdu]

        n_side = pars.n_side  # this is just for convenience since we'll use it so much

        for sca in range(1, 1 + pars.n_sca):
            im = np.zeros((n_side, n_side), dtype=np.uint16)
            valid = True
            try:
                with asdf.open(infile.format(sca)) as a:
                    im[:, :] = np.clip(np.rint(softbias + a["roman"]["data"] / dslope), 0, 2**16 - 1).astype(
                        np.uint16
                    )
                    start_time = a["roman"]["meta"]["exposure"]["start_time"]
                    mjd = a["roman"]["meta"]["ephemeris"]["time"]
            except FileNotFoundError:
                valid = False
            new_hdu = fits.ImageHDU(im)
            new_hdu.header["ISVALID"] = (valid, "Was this SCA found?")
            self.hdulist.append(new_hdu)

        # save collected metadata
        phdu.header["MJD"] = mjd
        phdu.header["TSTART"] = str(start_time)

        self.hdulist = fits.HDUList(hdulist) # save this as an HDUList

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
