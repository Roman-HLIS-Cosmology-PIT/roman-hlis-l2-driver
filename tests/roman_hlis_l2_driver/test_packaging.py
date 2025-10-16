import roman_hlis_l2_driver


def test_version():
    """Check to see that we can get the package version"""
    assert roman_hlis_l2_driver.__version__ is not None
