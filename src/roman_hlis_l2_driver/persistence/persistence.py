# persistence flagging driver script

import re
import json
import numpy as np

import persistence_flag as pf
import asdf as a
from astropy.io import fits as f

from pyimcom.config import Config
from pathlib import Path


def run(cfg: str, l2dir: str, delta_t_prime_max = 1200.0, signal_threshold = 20000.0):
  """
    Parameters
    ----------
    cfg: str
        File path to pyIMCOM config json
    l2dir: str
        #
    delta_t_prime_max: float
        #
    signal_threshold: float
        #

    Returns
    -------
    persistence_mask: np.array
        #
  """
  # open the config file
  with open(cfg, 'r') as file:
    cfg_file = json.load(file)
    print("loaded config file")

  
  obs_file_path = cfg_file["OBSFILE"]
  print(f"obs_file_path points to: {obs_file_path}")
  obs_file = f.open(obs_file_path)
  print("loaded fits file")
  date_list = obs_file[1].data["date"]
  exptime_list = obs_file[1].data["exptime"]
  
  file_directory = Path(cfg_file["INDATA"][0])
  print(f"file_directory points to: {file_directory}")

  for file in file_directory.iterdir():
    print(f"this file is: {file.name}")
    
    persistence_mask = np.full((4088,4088),False)
    
    search_format = '_(\d+)_(\d+)\.asdf$'
    matches = re.search(search_format,file.name)

    if matches is None:
      continue
    obsid = matches[1]
    sca = matches[2]

    print(f"found id: {obsid}, and sca: {sca}")

    delta_t = 0
    delta_t_prime = 0
    prev_ids = []

    
    while delta_t_prime <= delta_t_prime_max:
      prev_id, delta_t = pf.get_prev_obs(obs_file_path,id=int(obsid))
      # delta_t returned in days, convert to seconds
      delta_t *= 86400 # (24*60*60) to go from days to secs
      exptime = exptime_list[prev_id]
      delta_t_prime = delta_t + exptime
      print(f"delta_t: {delta_t}, delta_t_prime: {delta_t_prime}, exptime: {exptime}")
      l2_directory = Path(cfg_file["INDATA"][0][:-3]+"/")
      
      for l2_file in l2_directory.iterdir():
        search_format = '_(\d+)_(\d+)\.asdf$'
        l2_matches = re.search(search_format,l2_file.name)
        if l2_matches is None:
          continue
        if l2_matches[1] == prev_id:
          print(f"now looking in the L2 directory @ {l2_directory}")
          print(f"L2 file obsid {l2_matches[1]} matches current prev_id")
          print(f"looking at file: {l2_file.name}")
          current_file = a.open(l2_file)
          data_array = current_file['roman']['data']
          ny,nx = data_array.shape
  
          for j in range(ny):
            for i in range(nx):
              dn_per_s = data_array[j][i]
              dn = dn_per_s * exptime
              if dn >= signal_threshold:
                print(f'Pixel at (y={j},x={i}) has signal of {dn}DN')
                persistence_mask[j][i] = True

        else:
          continue

        
        print("breaking while loop by setting delta_t_prime to 1201")
        delta_t_prime = 1201
        
  return persistence_mask
