"""
Test functions for internal object finding.
"""

import numpy as np
from roman_hlis_l2_driver.starcatalogs.internal import encirc_center


def test_encirc_center():
    """Test function for encirc_center."""

    x_, y_ = np.meshgrid(np.arange(7), np.arange(8))
    use = (x_ - 3) ** 2 + (y_ - 3) ** 2 < 10.5
    use[3, 3:5] = False
    use[6, 5] = True
    use[0, 2] = False
    use[0, 4] = False
    print(use)

    # 0 values
    use2 = np.zeros_like(use)
    x, y, r = encirc_center(use2)
    assert abs(x + 1) < 1e-6
    assert abs(y + 1) < 1e-6
    assert abs(r + 1) < 1e-6
    print(x, y, r)

    # 1 value
    use2[5, 2] = True
    x, y, r = encirc_center(use2)
    assert abs(x - 2) < 1e-6
    assert abs(y - 5) < 1e-6
    assert abs(r) < 1e-6
    print(x, y, r)

    # 2 values
    use2[3, 3] = True
    x, y, r = encirc_center(use2)
    assert abs(x - 2.5) < 1e-6
    assert abs(y - 4) < 1e-6
    assert abs(r - 1.118033988749895) < 1e-6
    print(x, y, r)

    # General
    x, y, r = encirc_center(use)
    print(x, y, r)
    assert abs(x - 3.0454545454545454) < 1e-6
    assert abs(y - 3.3181818181818183) < 1e-6
    assert abs(r - 3.3184931360807237) < 1e-6
