import copy
import glob
import json
import os
import sys

import asdf
import numpy as np
from astropy.io import fits
from pyimcom import imdestripe
from pyimcom.config import Settings as Stn

from .destripe_setup import setup_all_files


def destripe_one_layer(cfg_file, noiseid=None, verbose=False):
    """
    Destripes one layer from the indicated set of files.

    The data in the files are *overwritten* and the ``processinfo`` leaf gets a new
    field: ``file["processinfo"]["destripe"]`` gives the number of
    noise layers that have been destriped, and ``file["processinfo"]["destripe_complete"]``
    is set to True if everything has been destriped.

    Parameters
    ----------
    cfg_file : str
        The configuration file.
    noiseid : int or None, optional
        If specified, destripes that particular noise layer.
        Otherwise, does the science layer.
    verbose : bool, optional
        Whether to talk a lot to the output.

    Returns
    -------
    int or None
        Number of noise layers. None if no files found.

    """

    # first get the file prefix and information
    with open(cfg_file, "r") as file:
        cfg = json.load(file)
    file_prefix = cfg["INDATA"][0]
    file_format = cfg["INDATA"][1]
    filter = Stn.RomanFilters[cfg["FILTER"]]
    out_prefix = cfg["DSOBSFILE"] + filter

    # add tails as needed
    if file_format == "L2_2506":
        file_prefix += f"/sim_L2_{filter:s}"
    else:
        print("Please add the new format to destripe.")
        raise ValueError(f"unsupported format in destripe: {file_format}")
    if verbose:
        print("File prefix =", file_prefix)

    use_files = setup_all_files(file_prefix, out_prefix, wcs_order=3, noiseid=noiseid, verbose=verbose)
    if verbose:
        print("Files selected:", use_files)

    # if there aren't any files, return None
    if len(use_files) == 0:
        return None

    # now we know there are some files
    with asdf.open(use_files[0][0][:-5] + "_noise.asdf") as a:
        n_noise_layer = np.shape(a["noise"])[0]
    if verbose:
        print("Number of noise layers:", n_noise_layer)

    # check we're in the right place
    if noiseid is not None:
        for fp in use_files:
            with asdf.open(fp[0], mode="rw") as a_in:
                if a_in["processinfo"]["destripe"] != noiseid:
                    raise ValueError(
                        f'Destriping counter not at the right noise field: '
                        f'expected {a_in["processinfo"]["destripe"]}, got {noiseid}'
                    )

    # cleanup output directory (except for overlap matrices)
    clearfiles = glob.glob(os.path.join(cfg["DSOUT"][0] + "/masks", "*_mask.fits"))
    clearfiles += glob.glob(os.path.join(cfg["DSOUT"][0] + "/masks", "*_mask.fits.lock"))
    clearfiles += glob.glob(os.path.join(cfg["DSOUT"][0], "*.fits"))
    if noiseid is None:
        # files to clear only the first time; these are reused for each noise layer
        clearfiles += glob.glob(os.path.join(cfg["DSOUT"][0], "ovmat.npy"))
        clearfiles += glob.glob(os.path.join(cfg["DSOUT"][0], "SCA_list.txt"))
        clearfiles += glob.glob(os.path.join(cfg["DSOUT"][0], "*.out"))
    if verbose:
        print("Clearing files:")
        for f in clearfiles:
            print("    ", f)
        print("")

    for p in clearfiles:
        os.remove(p)

    # main destriping
    dsout = imdestripe.main(cfg_file, overlaponly=False)
    if verbose:
        print("Output -->", dsout)
        print("")
        sys.stdout.flush()

    # now copy back
    for fp in use_files:
        if noiseid is None:
            with asdf.open(fp[0], mode="r", lazy_load=False) as a:
                a_in = copy.deepcopy(a.tree)
            a_in["processinfo"]["destripe"] = 0
            a_in["processinfo"]["destripe_complete"] = False
            with fits.open(dsout + fp[1] + ".fits") as f:
                a_in["destripe_orig"] = f[0].data.astype(np.float32)
        else:
            with asdf.open(fp[0], mode="r", lazy_load=False) as a:
                a_in = copy.deepcopy(a.tree)
            with fits.open(dsout + fp[1] + ".fits") as f:
                with asdf.open(fp[0][:-5] + "_noise.asdf", mode="rw") as anoise_in:
                    anoise_in_tree = copy.deepcopy(anoise_in.tree)
                anoise_in_tree["noise"][noiseid, :, :] = (f[0].data - a_in["destripe_orig"]).astype(
                    np.float16
                )
            asdf.AsdfFile(tree=anoise_in_tree).write_to(fp[0][:-5] + "_noise.asdf")
            a_in["processinfo"]["destripe"] += 1
            if a_in["processinfo"]["destripe"] == n_noise_layer:
                a_in["roman"]["data"] = np.copy(a_in["destripe_orig"])
                del a_in["destripe_orig"]
                a_in["processinfo"]["destripe_complete"] = True
        asdf.AsdfFile(tree=a_in).write_to(fp[0])

    return n_noise_layer


def destripe_all_layers(cfg_file, verbose=False):
    """
    Destripe all layer from the indicated set of files (including noise).

    The keyword ``file["processinfo"]["destripe_complete"]``
    is set to True if this executes correctly.

    Parameters
    ----------
    cfg_file : str
        The configuration file.
    verbose : bool, optional
        Whether to talk a lot to the output.

    Returns
    -------
    int
        Number of noise layers. None if no files found.

    """

    n_noise_layer = destripe_one_layer(cfg_file, verbose=verbose)
    for i_noise in range(n_noise_layer):
        destripe_one_layer(cfg_file, noiseid=i_noise, verbose=verbose)

    return n_noise_layer
