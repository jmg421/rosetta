"""Tests for rosetta.wrappers.phyloseq."""

import numpy as np
import pandas as pd
import pytest

from rosetta._errors import RDataError


def _phyloseq_available():
    try:
        from rosetta._deps import is_installed
        return is_installed("phyloseq")
    except Exception:
        return False


@pytest.fixture
def otu_table():
    np.random.seed(42)
    data = np.random.randint(0, 500, size=(5, 4))
    return pd.DataFrame(data, index=["OTU1", "OTU2", "OTU3", "OTU4", "OTU5"], columns=["S1", "S2", "S3", "S4"])


@pytest.fixture
def sample_meta():
    return pd.DataFrame({"site": ["gut", "gut", "skin", "skin"]}, index=["S1", "S2", "S3", "S4"])


def test_empty_otu_raises():
    from rosetta.wrappers.phyloseq import phyloseq
    with pytest.raises(RDataError, match="empty"):
        phyloseq(pd.DataFrame())


def test_negative_otu_raises(otu_table):
    from rosetta.wrappers.phyloseq import phyloseq
    bad = otu_table.copy()
    bad.iloc[0, 0] = -1
    with pytest.raises(RDataError, match="negative"):
        phyloseq(bad)


@pytest.mark.skipif(not _phyloseq_available(), reason="phyloseq not installed in R")
def test_phyloseq_creates_object(otu_table, sample_meta):
    from rosetta.wrappers.phyloseq import phyloseq
    ps = phyloseq(otu_table, sample_meta)
    assert ps is not None


@pytest.mark.skipif(not _phyloseq_available(), reason="phyloseq not installed in R")
def test_richness_returns_dataframe(otu_table, sample_meta):
    from rosetta.wrappers.phyloseq import phyloseq_richness
    result = phyloseq_richness(otu_table, sample_meta, measures=["Shannon", "Simpson"])
    assert isinstance(result, pd.DataFrame)
    assert "Shannon" in result.columns
    assert len(result) == 4  # 4 samples
