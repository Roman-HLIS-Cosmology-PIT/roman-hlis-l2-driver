# test_persistence_flag.py

import numpy as np
import pytest
from astropy.io import fits
from roman_hlis_l2_driver.persistence.persistence_flag import previous_obsid


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
    with pytest.raises(ValueError):
        result = previous_obsid(
            str(observation_table),
            row_number=2,
            obsid=62000.03000,
        )


def test_passing_no_search_argument_returns_none(observation_table):
    result = previous_obsid(str(observation_table))

    assert result is None
