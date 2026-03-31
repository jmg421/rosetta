"""DESeq2 differential expression wrapper."""

import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr

from .._bridge import _converter, to_r_dataframe, to_r_matrix, to_pandas, to_r_df
from .._deps import ensure_installed
from .._errors import RDataError, RFormulaError


def deseq2(counts: pd.DataFrame, metadata: pd.DataFrame, design: str = "~ condition", **kwargs) -> pd.DataFrame:
    """Run DESeq2 differential expression analysis.

    Args:
        counts: Gene count matrix (genes x samples) with non-negative integers.
        metadata: Sample metadata DataFrame with row names matching counts columns.
        design: R formula string for the experimental design.
        **kwargs: Additional arguments passed to DESeq2::DESeq().

    Returns:
        DataFrame with baseMean, log2FoldChange, lfcSE, stat, pvalue, padj.
    """
    if (counts < 0).any().any():
        raise RDataError("Count matrix contains negative values")
    if not set(counts.columns).issubset(set(metadata.index)):
        raise RDataError("Count matrix columns must match metadata row names")

    stats_pkg = importr("stats")
    with localconverter(_converter):
        try:
            stats_pkg.as_formula(design)
        except Exception as e:
            raise RFormulaError(f"Invalid design formula '{design}': {e}") from e

    ensure_installed("DESeq2")
    deseq2_pkg = importr("DESeq2")

    r_counts = to_r_matrix(counts)
    r_metadata = to_r_dataframe(metadata)
    r_design = ro.Formula(design)

    with localconverter(_converter):
        dds = deseq2_pkg.DESeqDataSetFromMatrix(
            countData=r_counts, colData=r_metadata, design=r_design,
        )
        dds = deseq2_pkg.DESeq(dds, **kwargs)
        res = deseq2_pkg.results(dds)

    return to_pandas(to_r_df(res))
