# script for flagging persistence,
# should be able to go back and look at the previous images from that SCA in order to flag.

import numpy as np
from astropy.io import fits as f

def previous_obsid(obsid: int, row_number: int, inputfile: str):
	""" Given an observation ID or row number in the table,
        returns the ID of the observation that came before it.
    Parameters
    ---------
    OBSID: int
        Observation ID
    row_number: int
        Row number of
    inputfile: str or str-like
        Path to observation table fits file
    Returns
    -------
    obsid_prev: int
        Observation ID of the previous observation to the one passed in.
	"""
	in_file = f.open(inputfile)
	
	dates = in_file[1].data["date"]
	sorted_dates = np.argsort(dates)
	
	current_obsid = np.where(sorted_dates == obsid)
	current_obsid_index = current_obsid[0][0]
	
	if current_obsid_index == 0:
		return None

	return dates[float(current_obsid_index) - 1]
