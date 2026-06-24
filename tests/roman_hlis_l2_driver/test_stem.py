"""Test for file name utilities."""

import pytest
from roman_hlis_l2_driver.name_util import stem_l2, stem_mask


def test_stem():
    """Test the stem_l2 function."""

    f = stem_l2("L2_2506", "J129")
    assert f == "/sim_L2_J129"

    with pytest.raises(ValueError):
        stem_l2("undefined_format", "J129")


def test_stem_mask():
    """Test the stem_mask function."""

    tail, hdu, bits = stem_mask("L2_2506")
    assert tail == "_mask.fits"
    assert hdu == "MASK"
    assert bits == 1

    with pytest.raises(ValueError):
        stem_mask("undefined_format")
