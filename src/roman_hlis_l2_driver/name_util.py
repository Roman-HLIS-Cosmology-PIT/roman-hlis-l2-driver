"""
Utilities for file names and formats.
"""


def stem_l2(file_format, filter):
    """
    Utility to return the stem format for an L2 file.

    Paramters
    ---------
    file_format : str
        The format type (same as used in PyIMCOM).
    filter : str
        The 4-character filter code.

    Returns
    -------
    str
        A string to append to a directory to get the stem.

    """

    if file_format == "L2_2506":
        return f"/sim_L2_{filter:s}"
    else:
        print("Please add the new format to name_util.")
        raise ValueError(f"unsupported format in name_util: {file_format}")
