# script for flagging persistence,
# should be able to go back and look at the previous images from that SCA in order to flag.

import numpy as np

def find_previous_sca(obsid = None, row_number = None, inputfile = None)
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


	return obsid_prev
