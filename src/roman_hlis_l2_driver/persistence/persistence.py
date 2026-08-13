# persistence flagging driver script

import re
import json
import numpy as np

import persistence_flag as pf
from astropy.io import fits as f
from pyimcom.config import Config


myfile = "/users/PCON0003/cond0007/Jul26/config-Jul26-H1.json"
cfg = Config(myfile)

def run(cfg, delta_t_prime_max = 1200.0, signal_threshold = 20000.0, l2dir: str):
  """
    Parameters
    ----------
    cfg: str
        File path to pyIMCOM config json
    delta_t_prime_max: float
        #
    signal_threshold: float
        #
    l2dir: str
        #
    Returns
    -------
    prev_id: int
        #
    """
  persistence_mask = np.full((4088,4088),False)
  
  with open(cfg, 'r') as file:
    cfg_file = json.load(file)

  obs_file_path = cfg_file["OBSFILE"]
  obs_file = f.open(obs_file_path)
  obs_list = obs_file[1].data["date"]
  exptime_list = obs_file[1].data["exptime"]
  
  file_directory = Path(cfg_file["INDATA"][0])

  for file in file_directory.iterdir():
    search_format = r'(\d+)_(\d+)\.asdf$'
    matches = re.search(search_format,cfg["OBSFILE"])
    obsid = matches[1]
    sca = matches[2]

    delta_t = 0
    delta_t_prime = 0
    prev_ids = []
    while delta_t_prime <= delta_t_prime_max:
      prev_id, delta_t = pf.get_prev_obs(obs_list,id=obsid)
      # delta_t returned in years, convert to seconds
      delta_t *= 86400 # (24*60*60) to go from years to secs
      delta_t_prime = delta_t + exptime[prev_id]
      prev_ids.append(prev_id)

  
  
  
