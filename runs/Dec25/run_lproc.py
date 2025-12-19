"""
Driver script for Level 2->2.1 processing.

Usage::

    python run_lproc.py <config> <asdf_mask>

Note that the asdf_mask is just used for visualization of the full set of
outliers masked, the mask to be passed to PyIMCOM is saved as an additional field
in the Level 2.1 files.

"""

import sys

# import numpy as np
# from astropy.io import fits
from roman_hlis_l2_driver.destripe_interface.destripe import destripe_all_layers
from roman_hlis_l2_driver.outliers.outlier_flagging import OutlierMap

# u = OutlierMap(sys.argv[1], max_workers=24)
# mask, out = u.outlier_mask(84)
# mask = mask.astype(np.int8)
# fits.PrimaryHDU(out).writeto("im1.fits", overwrite=True)
# fits.PrimaryHDU(mask).writeto("ma1.fits", overwrite=True)
# del u

destripe_all_layers(sys.argv[1], verbose=True)
OutlierMap(sys.argv[1], max_workers=27, run_and_save=sys.argv[2])
