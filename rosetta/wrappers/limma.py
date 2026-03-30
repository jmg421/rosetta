"""limma-voom differential expression wrapper."""

import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr

from .._bridge import _converter, to_r_matrix, to_r_dataframe, to_pandas
from .._deps import ensure_installed
from .._errors import RDataError, RFormulaError


def limma_voom(counts: pd.DataFrame, metadata: pd.DataFrame, design: str = "~ condition", **kwargs) -> pd.DataFrame:
    """Run limma-voom differential expression analysis.

    Args:
        counts: Gene count matrix (genes x samples) with non-negative integers.
        metadata: Sample metadata DataFrame with row names matching counts columns.
        design: R formula string for the experimental design.
        **kwargs: Additional arguments passed to lmFit().

    Returns:
        DataFrame with logFC, AveExpr, t, P.Value, adj.P.Val, B.
    """
    if (counts < 0).any().any():
        raise RDataError("Count matrix contains negative values")
    if not set(counts.columns).issubset(set(metadata.index)):
        raise RDataError("Count matrix columns must match metadata row names")

    with localconverter(_converter):
        try:
            ro.r["as.formula"](design)
        except Exception as e:
            raise RFormulaError(f"Invalid design formula '{design}': {e}") from e

    ensure_installed("limma")
    ensure_installed("edgeR")

    limma_pkg = importr("limma")
    edger_pkg = importr("edgeR")
    stats_pkg = importr("stats")

    r_counts = to_r_matrix(counts)
    r_metadata = to_r_dataframe(metadata)

    with localconverter(_converter):
        r_design_matrix = stats_pkg.model_matrix(ro.Formula(design), data=r_metadata)
        dge = edger_pkg.DGEList(counts=r_counts)
        dge = edger_pkg.calcNormFactors(dge)
        v = limma_pkg.voom(dge, r_design_matrix)
        fit = limma_pkg.lmFit(v, r_design_matrix, **kwargs)
        fit = limma_pkg.eBayes(fit)
        r_df = limma_pkg.topTable(fit, number=ro.r["nrow"](r_counts))

    return to_pandas(r_df)
