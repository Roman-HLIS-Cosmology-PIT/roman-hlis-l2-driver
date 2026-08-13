# persistence flagging driver script

import re
import json
import numpy as np

import persistence_flag as pf
import asdf as a
from astropy.io import fits as f

from pyimcom.config import Config
from pathlib import Path


def run(cfg: str, l2dir = None, delta_t_prime_max = 1200.0, signal_threshold = 20000.0):
    """
    Search previous L2 observations for pixels that may produce
    persistence in the current L2.1 observation.
    Parameters
    ----------
    cfg : str
        File path to pyIMCOM config JSON.
    l2dir : str
        Directory containing L2 ASDF files.
    delta_t_prime_max : float
        Maximum persistence lookback time in seconds.
    signal_threshold : float
        Minimum previous-observation signal in DN for a pixel
        to be considered capable of causing persistence.
    Returns
    -------
    mask_list : list
        List containing the cumulative persistence mask and the
        previous observation IDs used to construct it.
    """
    # open the config file
    with open(cfg, 'r') as file:
        cfg_file = json.load(file)
        
    print("loaded config file")

    obs_file_path = cfg_file["OBSFILE"]
    print(f"OBSFILE points to: {obs_file_path}")
    with f.open(obs_file_path) as obs_file:
        exptime_list = obs_file[1].data["exptime"].copy
        
    date_list = obs_file[1].data["date"]
    exptime_list = obs_file[1].data["exptime"]

    # L2.1 directory containing current observations
    file_directory = Path(cfg_file["INDATA"][0])
    print(f"INDATA: {file_directory}")

    # L2 directory containing previous observations
    if l2dir = None:
        l2_directory = Path(cfg_file["INDATA"][0][:-3]+"/")
    else:
        l2_directory = Path(l2dir)
    print(f"L2 directory: {l2_directory}")

    search_pattern = re.compile(r"_(\d+)_(\d+)\.asdf$")
    l2_files = {}

    for l2_file in l2_directory.iterdir():
        match = search_pattern.search(l2_file.name)
        if match is None:
            continue
        l2_obsid = int(match.group(1))
        l2_sca = int(match.group(2))
        l2_files[(l2_obsid, l2_sca)] = l2_file

    print(f"Indexed {len(l2_files)} L2 files")

    mask_list =[]

    for infile in file_directory.iterdir():
        matches = search_pattern.search(infile.name)
    if matches is None:
      continue
    current_obsid = int(matches.group(1))
    sca = int(matches.group(2))

    print()
    print(f"Current L2.1 file: {infile.name}")
    print(f"Current OBSID: {current_obsid}, SCA: {sca}")

    persistence_mask = None

    search_obsid = current_obsid
    total_delta_t = 0.0

    while True:
        print(f"Searching for observation before OBSID {search_obsid}")
        result = pf.get_prev_obs(obs_file_path,id=search_obsid)

        # Protect against reaching the beginning of the table.
        if result is None:
            print("No previous observation found.")
            break

        prev_id, delta_t = result

        prev_id = int(prev_id)
      
        # Convert time difference from days to secs
        delta_t_seconds = float(delta_t) * 86400.0

        # Since we are recursively stepping backwards,
        # accumulate the time differences.
        total_delta_t += delta_t_seconds
        exptime = float(exptime_list[prev_id])
        delta_t_prime = total_delta_t + exptime

        print(f"Previous OBSID: {prev_id}, Step delta_t:   {delta_t_seconds:.3f} sec")
        print(f"Exposure time:  {exptime:.3f} sec, delta_t_prime:  {delta_t_prime:.3f} sec")
        print(f"Total delta_t:  {total_delta_t:.3f} sec")

        if delta_t_prime > delta_t_prime_max:
            print(f"delta_t_prime = {delta_t_prime:.3f} sec > {delta_t_prime_max:.3f} sec")
            print("Persistence lookback window reached.")
            break

        prev_ids.append(prev_id)

        # Find matching L2 file
        key = (prev_id)

        if key not in l2_files:
            print(f"No L2 file found for OBSID {prev_id}, SCA {sca}")
            # Still move backwards to the next observation.
            search_obsid = prev_id
            continue

        l2_file = l2_files[key]
        print(f"Found matching L2 file: {l2_file}")
      
        # Read previous observation
        with a.open(l2_file) as current_file:
            data_array = np.asarray(current_file["roman"]["data"])
        # L2 is DN / second.
        signal_dn = data_array * exptime

        # Persistence mask for THIS previous observation
        previous_mask = signal_dn >= signal_threshold
        n_persistence_pixels = np.count_nonzero(previous_mask)

        print(f"{n_persistence_pixels} pixels exceed, {signal_threshold:.1f} DN")

        # Add this observation's pixels to cumulative mask
        if persistence_mask is None:
            persistence_mask = np.zeros_like(previous_mask,dtype=bool,)

        persistence_mask |= previous_mask
      
        # Continue recursive search
        search_obsid = prev_id

    # done searching for previous obs
    if persistence_mask is None:
        print("No qualifying previous observations found.")
        persistence_mask = np.zeros((4088, 4088),dtype=bool,)

    print(f"Previous observations searched: {prev_ids}")
    print(f"Total pixels in cumulative persistence mask: {np.count_nonzero(persistence_mask)}")

    print("DEBUG MODE: stopping after first valid L2.1 file.")
    return mask_list

return mask_list
