# script for flagging persistence,
# should be able to go back and look at the previous images from that SCA in order to flag.

import numpy as np
from astropy.io import fits as f


def previous_obsid(inputfile: str, row_number=None, obsid=None):
    """
        Given an observation ID (or a row number) and a path to a FITS file,
    returns the ID (or row number) of the observation that came before it
        in the file.

    Exactly one of `row_number` or `obsid` must be provided.

    Parameters
    ----------
    inputfile: str or str-like
        Path to observation table fits file.
    row_number: int, optional
        Unsorted row number, or index, of an observation within the fits file.
    obsid: int, optional
        Observation ID, taken as the modified Julian date.

        Returns
    -------
    int
            One of either "obsid_prev" or "prev_row", depending on whether `obsid` or
            `row_number` is provided:

                - ``obsid_prev`` : Observation ID of the previous observation to the one passed in.
        - ``prev_row`` : Row number of the observation previous to the row_number passed,
                were the rows sorted by date.

    """

    in_file = f.open(inputfile)
    all_obs = in_file[1].data["date"]

    if row_number is not None and obsid is not None:
        raise ValueError("Please pick EITHER an obsid or a row_number to search by, not both")

    if row_number is not None and obsid is None:
        sorted_index = int(np.where(np.argsort(all_obs) == row_number)[0][0])

        if sorted_index == 0:
            return None

        prev_row = int(np.argsort(all_obs)[sorted_index - 1])
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


def get_prev_obs(obstable: str, nbackup = 1, id = None):
    """ Find from obstable the id that is 'nbackup' previous
        Parameters
        ----------
        obstable: str
            File path for a FITS table of observations.
        nbackup: int
            Number of steps to take back (e.g., 1 is the previous observation).
        id: int
            The observation row # in the table.
        Returns
        -------
        prev_id: int
            The observation row # of the previous observation.
        delta_t: float
            Time elapsed from previous obs to this one (in days).
    """
    if id is None:
        raise ValueError("Please enter an id")
        
    in_file = f.open(obstable)
    all_obs = in_file[1].data["date"]

    current_obs = all_obs[id]
    delta_t = 0
    for backup_step in range(nbackup):
        prev_id = previous_obsid(inputfile = obstable, row_number = (id - backup_step) )

    prev_obs = all_obs[prev_id]
    delta_t = current_obs - prev_obs
    return prev_id, delta_t

def get_obs_date(obstable: str, id = None):
    """ Given an id, retrieve from obstable the date of the associated obs
        Parameters
        ----------
        obstable: str
            File path for a FITS table of observations.
        id: int
            The observation row # in the table.
        Returns
        -------
        obs_date: float
            Date of the requested obs, in modified Julian date.
    """
    in_file = f.open(obstable)
    all_obs = in_file[1].data["date"]

    obs_date = all_obs[id]
    return obs_date
