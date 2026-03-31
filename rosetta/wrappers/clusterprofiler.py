"""clusterProfiler gene set enrichment wrapper."""

from typing import List, Optional

import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr

from .._bridge import _converter, to_pandas, to_r_df
from .._deps import ensure_installed
from .._errors import RDataError


def enrichment(
    gene_list: List[str],
    organism: str = "org.Hs.eg.db",
    ont: str = "BP",
    pvalue_cutoff: float = 0.05,
    key_type: str = "SYMBOL",
    **kwargs,
) -> pd.DataFrame:
    """Run GO enrichment analysis via clusterProfiler::enrichGO.

    Args:
        gene_list: List of gene identifiers.
        organism: OrgDb annotation package (default: human).
        ont: GO ontology — "BP", "MF", or "CC".
        pvalue_cutoff: Adjusted p-value cutoff.
        key_type: Type of gene identifier (SYMBOL, ENTREZID, ENSEMBL, etc.).
        **kwargs: Additional arguments passed to clusterProfiler::enrichGO().

    Returns:
        DataFrame with enrichment results (ID, Description, pvalue, p.adjust, etc.).
    """
    if not gene_list:
        raise RDataError("gene_list must not be empty")

    ensure_installed("clusterProfiler")
    ensure_installed(organism)

    cp_pkg = importr("clusterProfiler")

    with localconverter(_converter):
        r_genes = ro.StrVector(gene_list)
        result = cp_pkg.enrichGO(
            gene=r_genes,
            OrgDb=organism,
            ont=ont,
            pvalueCutoff=pvalue_cutoff,
            keyType=key_type,
            **kwargs,
        )

    return to_pandas(to_r_df(result))
