"""Main driver function for imdestripe as called from Level 2-->2.1."""

import concurrent.futures
import copy
import glob
import json
import os
import sys

import asdf
import fitsio
import numpy as np
from pyimcom import imdestripe
from pyimcom.config import Settings as Stn

from ..name_util import stem_l2
from .destripe_setup import setup_all_files


def destripe_one_layer(cfg_file, noiseid=None, verbose=False, tracktest=False):
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
    tracktest : bool, optional
        If True, does one setup outside the parallelization for coverage tracking.
        You probably won't use this in general.

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

    use_files = setup_all_files(
        file_prefix, out_prefix, wcs_order=3, noiseid=noiseid, verbose=verbose, tracktest=tracktest
    )
    if verbose:
        print("Files selected:", use_files)

    # if there aren't any files, return None
    if len(use_files) == 0:
        return None

    # now we know there are some files
    with asdf.open(use_files[0][0][:-5] + "_noise.asdf", memmap=True) as a:
        n_noise_layer = np.shape(a["noise"])[0]
    if verbose:
        print("Number of noise layers:", n_noise_layer)
    nside = Stn.sca_nside

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
            # save the destriped image as a numpy memmap
            im = np.memmap(dsout + fp[1] + "_image.npy", dtype=np.float32, mode="w+", shape=(nside, nside))
            with fitsio.FITS(dsout + fp[1] + ".fits") as f:
                im[:, :] = f[0][:, :]
                im.flush()
        else:
            if noiseid == 0:
                noise_ds = np.memmap(
                    dsout + fp[1] + "_noise.npy",
                    dtype=np.float16,
                    mode="w+",
                    shape=(n_noise_layer, nside, nside),
                )
            else:
                noise_ds = np.memmap(
                    dsout + fp[1] + "_noise.npy",
                    dtype=np.float16,
                    mode="r+",
                    shape=(n_noise_layer, nside, nside),
                )
            im = np.memmap(dsout + fp[1] + "_image.npy", dtype=np.float32, mode="r", shape=(nside, nside))
            with fitsio.FITS(dsout + fp[1] + ".fits") as f:
                noise_ds[noiseid, :, :] = f[0][:, :] - im
                noise_ds.flush()
            if noiseid == n_noise_layer - 1:
                # last noise layer -- copy back, write as an ASDF file
                with asdf.open(fp[0], mode="r", lazy_load=False) as a:
                    a_in = copy.deepcopy(a.tree)
                a_in["roman"]["data"][:, :] = im
                a_in["processinfo"]["destripe_complete"] = True
                asdf.AsdfFile(tree=a_in).write_to(fp[0])
                with asdf.open(fp[0][:-5] + "_noise.asdf", mode="r", lazy_load=False) as a:
                    a_in = copy.deepcopy(a.tree)
                a_in["noise"][:, :, :] = noise_ds
                asdf.AsdfFile(tree=a_in).write_to(fp[0][:-5] + "_noise.asdf")
                # cleanup
                os.remove(dsout + fp[1] + "_noise.npy")
            del im, noise_ds

        # clean up image information
        if noiseid == n_noise_layer or n_noise_layer == 0:
            os.remove(dsout + fp[1] + "_image.npy")

    # execute one outside the executor for code coverage
    _cp(use_files[0])
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as e:
        iter = e.map(_cp, use_files[1:])
        while True:
            try:
                next(iter)
            except StopIteration:
                break  # End of iterator
            except Exception as exc:
                raise RuntimeError("file copyback error") from exc

    return n_noise_layer


def destripe_all_layers(cfg_file, verbose=False, tracktest=False):
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
    tracktest : bool, optional
        If True, does one setup outside the parallelization for coverage tracking.
        You probably won't use this in general.

    Returns
    -------
    int
        Number of noise layers. None if no files found.

    """

    n_noise_layer = destripe_one_layer(cfg_file, verbose=verbose, tracktest=tracktest)
    for i_noise in range(n_noise_layer):
        destripe_one_layer(cfg_file, noiseid=i_noise, verbose=verbose)

    return n_noise_layer
