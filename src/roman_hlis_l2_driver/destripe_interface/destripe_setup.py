import concurrent.futures
import os
import sys

import asdf
import numpy as np
from astropy.io import fits
from pyimcom.config import Settings as Stn
from pyimcom.wcsutil import LocWCS


def _setup_one_file(args):
    """
    Converts one file + ancillary information to FITS.
    """

    fp, outprefix, N, noiseid, wcs_order, verbose = args

    # convert the file to FITS, including 4th order TAN-SIP approximation
    outfile = outprefix + fp[1] + ".fits"
    with asdf.open(fp[0]) as input_file:
        # the WCS
        this_wcs = LocWCS(input_file["roman"]["meta"]["wcs"], N=N)  # Stn.sca_nside)
        this_wcs.wcs_approx_sip(p_order=wcs_order)

        # the main data
        if noiseid is None:
            phdu = fits.PrimaryHDU(
                input_file["roman"]["data"], header=this_wcs.approx_wcs.to_header(relax=True)
            )
        else:
            if verbose:
                print("Getting noise from", fp[0][:-5] + "_noise.asdf")
                sys.stdout.flush()
            with asdf.open(fp[0][:-5] + "_noise.asdf") as noise_file:
                phdu = fits.PrimaryHDU(
                    input_file["roman"]["data"] + noise_file["noise"][noiseid, :, :].astype(np.float32),
                    header=this_wcs.approx_wcs.to_header(relax=True),
                )

        # the mask
        mfname = fp[0][:-5] + "_mask.fits"
        with fits.open(mfname) as mf:
            if verbose:
                print(">>", outfile)
                print("    mask file:", mfname)
                print("    max WCS polynomial fit error:", this_wcs.wcs_max_err, "pix")
                sys.stdout.flush()
            fits.HDUList([phdu, mf[1]]).writeto(outfile, overwrite=True)


def setup_all_files(fprefix, outprefix, max_files=None, wcs_order=4, noiseid=None, verbose=False):
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
    noiseid : int, optional
        If specified, which noise layer to add from the
        fprefix + (numbers and underscores) + "_noise.asdf" file.
    verbose : bool, optional
        If True, print lots of data to the output.

    Returns
    -------
    list of (str, str)
        The list of file names selected and substrings ("obsid_sca").

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

    use_files2 = [(fp, outprefix, Stn.sca_nside, noiseid, wcs_order, verbose) for fp in use_files]
    with concurrent.futures.ProcessPoolExecutor() as e:
        e.map(_setup_one_file, use_files2)

    return use_files


# Simple case on OSC:
# python destripe_setup.py /fs/scratch/PCON0003/cond0007/fall2025/sim/L2/sim_L2_H158_ \
#     /fs/scratch/PCON0003/cond0007/temp/t_
if __name__ == "__main__":
    file_prefix = sys.argv[1]
    out_prefix = sys.argv[2]
    if len(sys.argv) > 3:
        setup_all_files(
            file_prefix, out_prefix, max_files=4, wcs_order=3, noiseid=int(sys.argv[3]), verbose=True
        )
    else:
        setup_all_files(file_prefix, out_prefix, max_files=4, wcs_order=3, verbose=True)
