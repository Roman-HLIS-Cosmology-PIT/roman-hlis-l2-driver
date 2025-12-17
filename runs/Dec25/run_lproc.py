"""
Driver script for Level 2->2.1 processing.

Usage::

    python run_lproc.py <config> <asdf_mask>

Note that the asdf_mask is just used for visualization of the full set of
outliers masked, the mask to be passed to PyIMCOM is saved as an additional field
in the Level 2.1 files.

"""

import sys

from roman_hlis_l2_driver.destripe_interface.destripe import destripe_all_layers
from roman_hlis_l2_driver.outliers.outlier_flagging import OutlierMap

destripe_all_layers(sys.argv[1], verbose=True)
OutlierMap(sys.argv[1], max_workers=12, run_and_save=sys.argv[2])
