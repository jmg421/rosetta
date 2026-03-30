"""Tests for rosetta.wrappers.clusterprofiler."""

import pandas as pd
import pytest

from rosetta._errors import RDataError


def _cp_available():
    try:
        from rosetta._deps import is_installed
        return is_installed("clusterProfiler") and is_installed("org.Hs.eg.db")
    except Exception:
        return False


def test_empty_gene_list_raises():
    from rosetta.wrappers.clusterprofiler import enrichment
    with pytest.raises(RDataError, match="empty"):
        enrichment([])


@pytest.mark.skipif(not _cp_available(), reason="clusterProfiler/org.Hs.eg.db not installed")
def test_enrichment_returns_dataframe():
    from rosetta.wrappers.clusterprofiler import enrichment
    # Well-known human genes likely to produce GO hits
    genes = ["TP53", "BRCA1", "EGFR", "MYC", "PTEN", "RB1", "AKT1", "KRAS", "BRAF", "PIK3CA"]
    result = enrichment(genes, pvalue_cutoff=0.5)
    assert isinstance(result, pd.DataFrame)
    if len(result) > 0:
        assert "Description" in result.columns
        assert "p.adjust" in result.columns
