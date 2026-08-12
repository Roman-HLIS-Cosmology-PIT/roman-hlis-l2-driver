# test_persistence_flag.py

import numpy as np
import pytest
from astropy.io import fits
from roman_hlis_l2_driver.persistence.persistence_flag import previous_obsid, get_prev_obs


@pytest.fixture
def observation_table(tmp_path):
    """
    Create a temporary FITS table whose observations are deliberately
    not stored in chronological order.

    Sort Order   Python index    DATE
    ----------   ------------    -----------
       2              0          62000.03000
       0              1          62000.01000
       3              2          62000.04000
       1              3          62000.02000
    Chronological order:
        row 1 -> row 3 -> row 0 -> row 2
    """
    dates = np.array([62000.03000, 62000.01000, 62000.04000, 62000.02000], dtype=np.float64)

    date_column = fits.Column(
        name="date",
        format="D",
        array=dates,
    )

    primary_hdu = fits.PrimaryHDU()
    table_hdu = fits.BinTableHDU.from_columns([date_column])

    filepath = tmp_path / "observation_table.fits"
    fits.HDUList([primary_hdu, table_hdu]).writeto(filepath)

    return filepath


def test_previous_obsid_exact_match(observation_table):
    """An exact OBSID match should return the next-lowest OBSID."""
    result = previous_obsid(
        str(observation_table),
        obsid=62000.03000,
    )

    assert result == pytest.approx(62000.02000)


def test_previous_obsid_between_observations(observation_table):
    """
    An OBSID between two recorded observations should return the
    largest recorded OBSID below it.
    """
    result = previous_obsid(
        str(observation_table),
        obsid=62000.03500,
    )

    assert result == pytest.approx(62000.03000)


def test_previous_obsid_after_last_observation(observation_table):
    """An OBSID after every observation should return the final OBSID."""
    result = previous_obsid(
        str(observation_table),
        obsid=62000.05000,
    )

    assert result == pytest.approx(62000.04000)


@pytest.mark.parametrize(
    "obsid",
    [
        62000.01000,  # Exactly the earliest observation
        62000.00000,  # Earlier than every observation
    ],
)
def test_previous_obsid_returns_none_when_no_previous_exists(
    observation_table,
    obsid,
):
    """Checks that None is passed if the observation_table does
    not contain an obsid previous to the one passed.
    """
    result = previous_obsid(
        str(observation_table),
        obsid=obsid,
    )

    assert result is None


def test_previous_row_number(observation_table):
    """
    As passed, row_number = 3 corresponds to DATE=62000.02000.

    The previous chronological observation is DATE=62000.01000,
    which is located in FITS row 1.
    """
    result = previous_obsid(
        str(observation_table),
        row_number=3,
    )

    assert result == 1


def test_first_chronological_row_has_no_previous_row(observation_table):
    """
    FITS row 1 contains the earliest observation, so it has no
    chronologically previous observation.
    """
    result = previous_obsid(
        str(observation_table),
        row_number=1,
    )

    assert result is None


def test_passing_both_search_arguments_prints_warning(
    observation_table,
):
    """Checks that ValueError is raised if we provide both arguments."""
    with pytest.raises(ValueError):
        previous_obsid(
            str(observation_table),
            row_number=2,
            obsid=62000.03000,
        )


def test_passing_no_search_argument_returns_none(observation_table):
    """Checks to make sure None is passed when neither
    row_number or obsid is passed
    """
    result = previous_obsid(str(observation_table))

    assert result is None


def test_get_prev_obs_returns_none_when_no_previous_exists(observation_table):
    """Checks that Nones are passed if the observation_table does
    not contain an id previous to the one passed.
    """
    result1, result2 = get_prev_obs(str(observation_table),id=1)

    assert result1 is None and result2 is None


def test_get_prev_obs_nbackup_not_supplied(observation_table):
    """
    As passed, row_number = 2 corresponds to DATE=62000.04000.

    The previous chronological observation is DATE=62000.03000,
    which is located in FITS row 0, and their difference is 0.01.
    """
    result1, result2 = get_prev_obs(str(observation_table), id=2)

    assert result1 == 0 and result2 == np.abs(result2 - 0.01) < 1.0e-5


def test_get_prev_obs_with_nbackup_passed(observation_table):
    """
    As passed, row_number = 0 corresponds to DATE=62000.03000.

    The observation that is (n=2) before is DATE=62000.01000,
    which is located in FITS row 1, and their difference is 0.02.
    """
    result1, result2 = get_prev_obs(str(observation_table), id=0, nbackup=2)

    assert result1 == 1 and result2 == np.abs(result2 - 0.02) < 1.0e-5
