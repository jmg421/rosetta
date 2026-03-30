"""Shared test fixtures."""

import numpy as np
import pandas as pd
import pytest


@pytest.fixture
def sample_counts():
    """3 genes x 4 samples count matrix."""
    np.random.seed(42)
    data = np.random.randint(0, 1000, size=(3, 4))
    return pd.DataFrame(data, index=["GeneA", "GeneB", "GeneC"], columns=["S1", "S2", "S3", "S4"])


@pytest.fixture
def sample_metadata():
    """Sample metadata with condition column."""
    return pd.DataFrame({"condition": ["control", "control", "treated", "treated"]}, index=["S1", "S2", "S3", "S4"])
