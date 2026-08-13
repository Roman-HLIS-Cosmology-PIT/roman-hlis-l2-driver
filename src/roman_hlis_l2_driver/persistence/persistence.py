# persistence flagging driver script

import re
import json
import numpy as np

import persistence_flag as pf
from astropy.io import fits as f
from pyimcom.config import Config


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
    prev_id: int
        #
    """
  # open the config file
  with open(cfg, 'r') as file:
    cfg_file = json.load(file)
    print("loaded config file")

  
  obs_file_path = cfg_file["OBSFILE"]
  obs_file = f.open(obs_file_path)
  print("loaded fits file")
  date_list = obs_file[1].data["date"]
  exptime_list = obs_file[1].data["exptime"]
  
  file_directory = Path(cfg_file["INDATA"][0])
  print(f"file_directory points to: {file_directory}")

  for file in file_directory.iterdir():
    persistence_mask = np.full((4088,4088),False)
    
    search_format = r'(\d+)_(\d+)\.asdf$'
    matches = re.search(search_format,file)
    obsid = matches[1]
    sca = matches[2]

    print(f"this file is: {file}, id: {obsid}, sca: {sca}")

    delta_t = 0
    delta_t_prime = 0
    prev_ids = []

    break

  
    
    while delta_t_prime <= delta_t_prime_max:
      prev_id, delta_t = pf.get_prev_obs(date_list,id=obsid)
      # delta_t returned in days, convert to seconds
      delta_t *= 86400 # (24*60*60) to go from days to secs
      delta_t_prime = delta_t + exptime[prev_id]
      l2_directory = Path(cfg_file["INDATA"][0][:-3]+"/")
      for file in l2_directory.iterdir():

  
  
  
