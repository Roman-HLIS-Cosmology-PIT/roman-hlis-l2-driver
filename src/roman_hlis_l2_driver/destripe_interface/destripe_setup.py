import os
import sys

import asdf
from astropy.io import fits
from pyimcom.config import Settings as Stn
from pyimcom.wcsutil import LocWCS


def setup_all_files(fprefix, outprefix, max_files=None, wcs_order=4, verbose=False):
    """
    Gets all the files starting with the specified format.

    The outputs are written as fits files starting with `outprefix`.

    Parameters
    ----------
    fprefix : str
        The prefix for the file names. Files will be accepted if they are of the
        format fprefix + (numbers and underscores) + ".asdf".
    outprefix : str
        The prefix for the output file names. Files will be written to
        outprefix + (same numbers and undercores) + ".fits".
    max_files : int or None, optional
        If provided, sets a maximum number of files (for testing).
    wcs_order : int, optional
        The polynomial fitting order for the TAN-SIP WCS.
    verbose : bool, optional
        If True, print lots of data to the output.

    Returns
    -------
    None

    """

    # get the input files to use. here use_files is a list of ordered pairs
    # (file, identifier)
    fdir, fileprefix = os.path.split(fprefix)
    n = len(fileprefix)
    numus = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "_"}
    use_files = []
    for f in os.listdir(fdir):
        if len(use_files) == max_files:
            break
        if f[:n] == fileprefix and f[-5:] == ".asdf" and all(c in numus for c in f[n:-5]):
            use_files.append((os.path.join(fdir, f), f[n:-5]))

    for fp in use_files:
        # convert the file to FITS, including 4th order TAN-SIP approximation
        outfile = outprefix + fp[1] + ".fits"
        with asdf.open(fp[0]) as input_file:
            # the WCS
            this_wcs = LocWCS(input_file["roman"]["meta"]["wcs"], N=Stn.sca_nside)
            this_wcs.wcs_approx_sip(p_order=wcs_order)

            # the main data
            phdu = fits.PrimaryHDU(
                input_file["roman"]["data"], header=this_wcs.approx_wcs.to_header(relax=True)
            )

            # the mask
            mfname = fp[0][:-5] + "_mask.fits"
            with fits.open(mfname) as mf:
                if verbose:
                    print(">>", outfile)
                    print("    mask file:", mfname)
                    print("    max WCS polynomial fit error:", this_wcs.wcs_max_err, "pix")
                fits.HDUList([phdu, mf[1]]).writeto(outfile, overwrite=True)


# Simple case on OSC:
# python destripe_setup.py /fs/scratch/PCON0003/cond0007/fall2025/sim/L2/sim_L2_H158_ \
#     /fs/scratch/PCON0003/cond0007/temp/t_
if __name__ == "__main__":
    file_prefix = sys.argv[1]
    out_prefix = sys.argv[2]
    setup_all_files(file_prefix, out_prefix, max_files=4, wcs_order=3, verbose=True)
