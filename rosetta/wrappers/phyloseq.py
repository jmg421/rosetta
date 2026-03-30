"""phyloseq microbiome analysis wrapper."""

from typing import Optional

import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr

from .._bridge import _converter, to_r_matrix, to_r_dataframe, to_pandas
from .._deps import ensure_installed
from .._errors import RDataError


def phyloseq(
    otu_table: pd.DataFrame,
    sample_data: Optional[pd.DataFrame] = None,
    tax_table: Optional[pd.DataFrame] = None,
) -> object:
    """Create a phyloseq object from pandas DataFrames.

    Args:
        otu_table: OTU/ASV count matrix (taxa x samples).
        sample_data: Sample metadata DataFrame.
        tax_table: Taxonomy table DataFrame (taxa x ranks).

    Returns:
        rpy2 phyloseq object for further analysis.
    """
    if otu_table.empty:
        raise RDataError("otu_table must not be empty")
    if (otu_table < 0).any().any():
        raise RDataError("OTU table contains negative values")

    ensure_installed("phyloseq")
    ps_pkg = importr("phyloseq")

    r_otu = to_r_matrix(otu_table)

    with localconverter(_converter):
        ps_otu = ps_pkg.otu_table(r_otu, taxa_are_rows=True)
        components = [ps_otu]

        if sample_data is not None:
            r_sd = to_r_dataframe(sample_data)
            ps_sd = ps_pkg.sample_data(r_sd)
            components.append(ps_sd)

        if tax_table is not None:
            r_tax = to_r_matrix(tax_table)
            ps_tax = ps_pkg.tax_table(r_tax)
            components.append(ps_tax)

        ps_obj = ps_pkg.phyloseq(*components)

    return ps_obj


def phyloseq_richness(
    otu_table: pd.DataFrame,
    sample_data: Optional[pd.DataFrame] = None,
    measures: Optional[list] = None,
) -> pd.DataFrame:
    """Calculate alpha diversity metrics.

    Args:
        otu_table: OTU/ASV count matrix (taxa x samples).
        sample_data: Sample metadata DataFrame.
        measures: Diversity measures (e.g. ["Shannon", "Simpson", "Chao1"]).

    Returns:
        DataFrame with diversity metrics per sample.
    """
    ensure_installed("phyloseq")
    ps_pkg = importr("phyloseq")

    ps_obj = phyloseq(otu_table, sample_data)

    with localconverter(_converter):
        kwargs = {}
        if measures:
            kwargs["measures"] = ro.StrVector(measures)
        r_df = ps_pkg.estimate_richness(ps_obj, **kwargs)

    return to_pandas(r_df)
