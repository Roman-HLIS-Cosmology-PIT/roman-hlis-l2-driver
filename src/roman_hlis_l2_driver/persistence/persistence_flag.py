# script for flagging persistence,
# should be able to go back and look at the previous images from that SCA in order to flag.

import numpy as np
from astropy.io import fits as f

def previous_obsid(obsid: int, inputfile: str):
	""" Given an observation ID or row number in the table,
        returns the ID of the observation that came before it.
    Parameters
    ---------
    obsid: int
        Observation ID
    inputfile: str or str-like
        Path to observation table fits file
    Returns
    -------
    obsid_prev: int
        Observation ID of the previous observation to the one passed in.
	"""
	in_file = f.open(inputfile)
	
	all_obsid = in_file[1].data["date"]
	sorted_obsid_indices = np.argsort(all_obsid)
	
	previous_obsid = None

	for current_index in sorted_obsid_indices:
		current_obsid = all_obsid[current_index]
		if current_obsid >= obsid:
			break
		previous_obsid = float(current_obsid)

	return previous_obsid
