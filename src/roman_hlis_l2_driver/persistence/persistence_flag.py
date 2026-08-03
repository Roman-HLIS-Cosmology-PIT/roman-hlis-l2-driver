# script for flagging persistence,
# should be able to go back and look at the previous images from that SCA in order to flag.

import numpy as np
from astropy.io import fits as f

def previous_obsid(inputfile: str, row_number = None, obsid = None):
	""" Given an observation ID (or a row number) and a path to a FITS file,
        returns the ID (or row number) of the observation that came before it
		in the file.
    Parameters
    ---------
    inputfile: str or str-like
        Path to observation table fits file.
    row_number: int
        Unsorted row number, or index, of an observation within the fits file.
    obsid: int
        Observation ID, taken as the modified Julian date.
    Returns
    -------
    obsid_prev: int
        Observation ID of the previous observation to the one passed in.
    prev_row: int
        Row number of the observation previous to the row_number passed,
		sorted chronologically.
	"""
	in_file = f.open(inputfile)
	all_obs = in_file[1].data["date"]

	if row_number is not None and obsid is not None:
		print("Please pick EITHER an obsid or a row_number to search by, not both")
	
	if row_number is not None and obsid is None:
		selected_row = int(np.where(np.argsort(all_obs) == row_number)[0][0])
		prev_row = int(np.argsort(all_obs)[selected_row - 1])
		
		if prev_row == 0:
			return None
			
		return prev_row

	if row_number is None and obsid is not None: 
		sorted_obsid_indices = np.argsort(all_obs)
		
		previous_obsid = None
	
		for current_index in sorted_obsid_indices:
			current_obsid = all_obs[current_index]
			if current_obsid >= obsid:
				break
			previous_obsid = float(current_obsid)
	
		return previous_obsid
