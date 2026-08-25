# persistence flagging driver script

import json
import re
import warnings
from pathlib import Path

import asdf as a
import numpy as np
from astropy.io import fits as f
from erfa import ErfaWarning

from . import persistence_flag as pf


def run(cfg: str, l2dir=None, delta_t_prime_max=1200.0, signal_threshold=20000.0, nmax=20):
    """
    Search previous L2 observations for pixels that may produce
    persistence in the current L2.1 observation.

    Parameters
    ----------
    cfg : str
        File path to pyIMCOM config JSON.
    l2dir : str, optional
        Directory containing L2 ASDF files.
    delta_t_prime_max : float
        Maximum persistence lookback time in seconds.
    signal_threshold : float
        Minimum previous-observation signal in DN for a pixel
        to be considered capable of causing persistence.
    nmax : int
        Maximum number of steps allowed to take back from L2.1 obsid.

    Returns
    -------
    mask_list : list
        List containing tuples of:
        - the cumulative persistence mask: np.array of boolean values
        - list of previous obsids: list of int
        - the obsid of the associated L2.1 file: int
    """
    # Load config file
    with open(cfg, "r") as file:
        cfg_file = json.load(file)

    print("loaded config file")

    obs_file_path = cfg_file["OBSFILE"]
    print(f"OBSFILE points to: {obs_file_path}")

    # Read what we need while the FITS file is open
    with f.open(obs_file_path) as obs_file:
        exptime_list = obs_file[1].data["exptime"].copy()

    # L2.1 directory containing current observations
    file_directory = Path(cfg_file["INDATA"][0])
    print(f"INDATA: {file_directory}")

    # L2 directory containing previous observations
    l2_directory = Path(cfg_file["INDATA"][0][:-3] + "/") if l2dir is None else Path(l2dir)

    print(f"L2 directory: {l2_directory}")

    # File-name pattern
    search_pattern = re.compile(r"_(\d+)_(\d+)\.asdf$")

    l2_files = {}
    for l2_file in l2_directory.iterdir():
        match = search_pattern.search(l2_file.name)
        if match is None:
            continue
        l2_obsid = int(match.group(1))

        # Save the first file encountered for this OBSID
        if l2_obsid not in l2_files:
            l2_files[l2_obsid] = l2_file

    print(f"Indexed {len(l2_files)} L2 observation IDs")

    mask_list = []

    # loop over L2.1 files
    for infile in sorted(file_directory.iterdir()):
        matches = search_pattern.search(infile.name)
        if matches is None:
            continue

        current_obsid = int(matches.group(1))
        sca = int(matches.group(2))

        print(f"Current L2.1 file: {infile.name}")
        print(f"Current OBSID: {current_obsid}, SCA: {sca}")

        persistence_mask = None

        # Keep track of the previous observations we actually consider
        prev_ids = []

        search_obsid = current_obsid
        total_delta_t = 0.0
        nmax_counter = 0

        # walk backwards through obsids
        while True:
            nmax_counter += 1
            print(f"Searching for observation before OBSID {search_obsid}")

            result = pf.get_prev_obs(obs_file_path, id=search_obsid)

            prev_id, delta_t = result
            # get_prev_obs() returns (None, None) when there is no
            # no previous observation to be found.
            if prev_id is None:
                print("No previous observation found.")
                break

            prev_id = int(prev_id)

            # Convert days -> seconds
            delta_t_seconds = float(delta_t) * 86400.0

            # Accumulate time as we step backward
            total_delta_t += delta_t_seconds
            exptime = float(exptime_list[prev_id])
            delta_t_prime = total_delta_t + exptime

            print(f"Previous OBSID: {prev_id}, Step delta_t: {delta_t_seconds:.3f} sec")
            print(f"Exposure time: {exptime:.3f} sec, delta_t_prime: {delta_t_prime:.3f} sec")
            print(f"Total delta_t: {total_delta_t:.3f} sec")

            # Persistence time-window test
            if delta_t_prime > delta_t_prime_max:
                print(f"delta_t_prime = {delta_t_prime:.3f} sec > {delta_t_prime_max:.3f} sec")
                print("Persistence lookback window reached.")
                break
            if nmax_counter > nmax:
                warnings.warn("Default nmax counter exceeding, exiting early.")
                break
            prev_ids.append(prev_id)

            # Find matching L2 file using OBSID only
            if prev_id not in l2_files:
                print(f"No L2 file found for OBSID {prev_id}")
                # Keep walking backward anyway
                search_obsid = prev_id
                continue

            l2_file = l2_files[prev_id]
            print(f"Found matching L2 file: {l2_file}")

            # Read previous observation
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    "ignore", message=r'ERFA function "dtf2d".*dubious year.*', category=ErfaWarning
                )
                with a.open(l2_file) as current_file:
                    data_array = np.asarray(current_file["roman"]["data"])

            # L2 data are DN / second
            signal_dn = data_array * exptime

            # Pixels capable of causing persistence
            previous_mask = signal_dn >= signal_threshold

            n_persistence_pixels = np.count_nonzero(previous_mask)

            print(f"{n_persistence_pixels} pixels exceed {signal_threshold:.1f} DN")

            # Add to cumulative persistence mask
            if persistence_mask is None:
                persistence_mask = np.zeros_like(previous_mask, dtype=bool)

            persistence_mask |= previous_mask

            # Move backward another observation
            search_obsid = prev_id

        if persistence_mask is None:
            print("No qualifying previous observations found.")
            print(f"Previous observations searched: {prev_ids}")
            continue

        print(f"Previous observations searched: {prev_ids}")
        print(f"Total pixels in cumulative persistence mask: {np.count_nonzero(persistence_mask)}")

        # Actually put our result into mask_list
        mask_list.append((persistence_mask, prev_ids, current_obsid))

    return mask_list
