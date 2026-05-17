"""Test for file name utilities."""

import pytest

from roman_hlis_l2_driver.name_util import stem_l2

def test_stem():
    """Test the stem_l2 function."""

    f = stem_l2("L2_2506", "J129")
    assert f == "/sim_L2_J129"

    with pytest.raises(ValueError):
        g = stem_l2("undefined_format", "J129")
