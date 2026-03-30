"""Tests for rosetta._bridge."""

import numpy as np
import pandas as pd
import pytest

from rosetta._bridge import to_r_matrix, to_r_dataframe, to_pandas
from rosetta._errors import RDataError


def test_dataframe_roundtrip(sample_counts):
    r_df = to_r_dataframe(sample_counts)
    result = to_pandas(r_df)
    pd.testing.assert_frame_equal(result, sample_counts.astype(float), check_dtype=False)


def test_matrix_roundtrip(sample_counts):
    r_mat = to_r_matrix(sample_counts)
    result = to_pandas(r_mat)
    np.testing.assert_array_equal(result, sample_counts.values)


def test_to_r_matrix_rejects_non_dataframe():
    with pytest.raises(RDataError):
        to_r_matrix("not a dataframe")


def test_to_r_dataframe_rejects_non_dataframe():
    with pytest.raises(RDataError):
        to_r_dataframe([1, 2, 3])
