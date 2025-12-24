"""Tests for some of the utility functions."""

import numpy as np
from roman_hlis_l2_driver.fullfov import get_t_eff

def test_get_t_eff():
    """Simple test function for timing."""

    K = np.array([-1, 0, 1]).astype(np.float32)/6.
    tau = 3.0 * np.arange(2,5).astype(np.float32)
    tbar = 3.0 * np.arange(2,5).astype(np.float32)

    t = get_t_eff(K, tau, tbar)
    assert np.abs(t-6.)<1e-3
