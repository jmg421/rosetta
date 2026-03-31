"""Seurat single-cell analysis wrapper."""

import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr

from .._bridge import _converter, to_r_matrix, to_pandas, to_r_df
from .._deps import ensure_installed
from .._errors import RDataError


def seurat(
    counts: pd.DataFrame,
    min_cells: int = 3,
    min_features: int = 200,
    n_variable_features: int = 2000,
    n_pcs: int = 10,
    resolution: float = 0.5,
    **kwargs,
) -> dict:
    """Run standard Seurat single-cell clustering pipeline.

    Args:
        counts: Raw count matrix (genes x cells).
        min_cells: Minimum cells a gene must appear in.
        min_features: Minimum genes a cell must express.
        n_variable_features: Number of variable features to select.
        n_pcs: Number of PCs for neighbors/clustering.
        resolution: Clustering resolution.
        **kwargs: Additional arguments passed to Seurat::FindClusters().

    Returns:
        Dict with:
            - "clusters": Series mapping cell barcodes to cluster IDs
            - "umap": DataFrame with UMAP coordinates
            - "variable_features": list of variable gene names
    """
    if counts.empty:
        raise RDataError("counts must not be empty")
    if (counts < 0).any().any():
        raise RDataError("Count matrix contains negative values")

    ensure_installed("Seurat")

    seurat_pkg = importr("Seurat")
    sobj_pkg = importr("SeuratObject")
    base_pkg = importr("base")
    methods_pkg = importr("methods")
    r_counts = to_r_matrix(counts)
    dims = ro.IntVector(range(1, n_pcs + 1))

    with localconverter(_converter):
        obj = sobj_pkg.CreateSeuratObject(
            counts=r_counts, **{"min.cells": min_cells, "min.features": min_features},
        )
        obj = seurat_pkg.NormalizeData(obj, verbose=False)
        obj = seurat_pkg.FindVariableFeatures(obj, nfeatures=n_variable_features, verbose=False)
        obj = seurat_pkg.ScaleData(obj, verbose=False)
        obj = seurat_pkg.RunPCA(obj, npcs=n_pcs, verbose=False)
        obj = seurat_pkg.FindNeighbors(obj, dims=dims, verbose=False)
        obj = seurat_pkg.FindClusters(obj, resolution=resolution, verbose=False, **kwargs)
        obj = seurat_pkg.RunUMAP(obj, dims=dims, verbose=False)

        meta = to_r_df(methods_pkg.slot(obj, "meta.data"))
        meta_df = to_pandas(meta)
        clusters = meta_df["seurat_clusters"]

        embeddings = to_r_df(sobj_pkg.Embeddings(obj, reduction="umap"))
        umap_df = to_pandas(embeddings)

        var_features = list(sobj_pkg.VariableFeatures(obj))

    return {
        "clusters": clusters,
        "umap": umap_df,
        "variable_features": var_features,
    }
