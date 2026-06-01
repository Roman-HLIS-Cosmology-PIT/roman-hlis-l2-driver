"""Main driver function for imdestripe as called from Level 2-->2.1."""

import concurrent.futures
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

from ..name_util import stem_l2
from .destripe_setup import setup_all_files


def destripe_one_layer(cfg_file, noiseid=None, verbose=False):
    """
    Destripes one layer from the indicated set of files.

    The data in the files are *overwritten* and the ``processinfo`` leaf gets a new
    field: ``file["processinfo"]["destripe_complete"]``
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

    Notes
    -----
    A previous version of this function placed a progress indicator in the ASDF file,
    but we turned this off because it was taking too long.

    """

    # first get the file prefix and information
    with open(cfg_file, "r") as file:
        cfg = json.load(file)
    file_prefix = cfg["INDATA"][0]
    file_format = cfg["INDATA"][1]
    filter = Stn.RomanFilters[cfg["FILTER"]]
    out_prefix = cfg["DSOBSFILE"] + filter

    # add tails as needed
    file_prefix += stem_l2(file_format, filter)
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
    def _cp(fp):
        if verbose:
            print("copy back", fp)
            sys.stdout.flush()
        if noiseid is None:
            with asdf.open(fp[0], mode="r", lazy_load=False) as a:
                a_in = copy.deepcopy(a.tree)
            a_in["processinfo"]["destripe"] = 0
            a_in["processinfo"]["destripe_complete"] = False
            with fits.open(dsout + fp[1] + ".fits") as f:
                a_in["destripe_orig"] = f[0].data.astype(np.float32)
            asdf.AsdfFile(tree=a_in).write_to(fp[0])
        else:
            with fits.open(dsout + fp[1] + ".fits") as f:
                with asdf.open(fp[0][:-5] + "_noise.asdf", mode="rw", memmap=True) as anoise_in:
                    with asdf.open(fp[0], memmap=True) as aorig:
                        anoise_in["noise"][noiseid, :, :] = (f[0].data - aorig["destripe_orig"]).astype(
                            np.float16
                        )
                    anoise_in.update()
            if noiseid == n_noise_layer - 1:
                # last noise layer
                with asdf.open(fp[0], mode="r", lazy_load=False) as a:
                    a_in = copy.deepcopy(a.tree)
                a_in["roman"]["data"] = np.copy(a_in["destripe_orig"])
                del a_in["destripe_orig"]
                a_in["processinfo"]["destripe_complete"] = True
                asdf.AsdfFile(tree=a_in).write_to(fp[0])
        return 0

    # execute one outside the executor for code coverage
    _cp(use_files[0])
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as e:
        e.map(use_files[1:])

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
