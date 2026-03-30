# rosetta — Technical Specification

## Overview

`rosetta` is a Python library that wraps R bioinformatics packages (DESeq2, edgeR, limma, etc.) via `rpy2`, providing a pandas-native API.

## Architecture

```
rosetta/
├── __init__.py          # Public API (rb.deseq2, rb.edger, etc.)
├── _bridge.py           # rpy2 session management and R↔Python type conversion
├── _deps.py             # R/Bioconductor package detection and installation
├── _errors.py           # R error translation to Python exceptions
├── wrappers/
│   ├── __init__.py
│   ├── deseq2.py        # DESeq2 wrapper
│   ├── edger.py         # edgeR wrapper
│   ├── limma.py         # limma-voom wrapper
│   ├── seurat.py        # Seurat wrapper
│   ├── clusterprofiler.py
│   └── phyloseq.py
└── tests/
    ├── conftest.py      # Shared fixtures (sample counts, metadata)
    ├── test_bridge.py
    └── test_deseq2.py
```

## Core Components

### 1. Bridge Layer (`_bridge.py`)

Manages a single `rpy2` R session and handles all type conversion.

Key conversions:
- `pandas.DataFrame` ↔ `R data.frame`
- `numpy.ndarray` ↔ `R matrix`
- Python `dict` ↔ `R named list`
- `None` ↔ `R NULL`

### 2. Dependency Manager (`_deps.py`)

On first use of a wrapper, checks if the required R package is installed. If missing, installs via `BiocManager::install()` with user confirmation.

### 3. Error Translation (`_errors.py`)

Catches `rpy2.rinterface_lib.embedded.RRuntimeError` and maps common R errors to descriptive Python exceptions:
- `RPackageMissing` — R package not installed
- `RFormulaError` — invalid design formula
- `RDataError` — incompatible input data (e.g. negative counts for DESeq2)

### 4. Wrapper Pattern

Each wrapper follows the same structure:

```python
def deseq2(counts: pd.DataFrame, metadata: pd.DataFrame, design: str, **kwargs) -> pd.DataFrame:
    """Run DESeq2 differential expression analysis."""
    ensure_installed("DESeq2")
    r_counts = to_r_matrix(counts)
    r_metadata = to_r_dataframe(metadata)
    # Call R functions via rpy2
    # Return results as pandas DataFrame
```

All wrappers:
- Accept pandas DataFrames as input
- Return pandas DataFrames as output
- Expose R parameters as Python keyword arguments
- Validate inputs before crossing the R boundary

## Design Decisions

- **Wrap, don't reimplement** — statistical correctness comes from the original R packages
- **Lazy R initialization** — R session starts on first wrapper call, not on import
- **One function per analysis** — `rb.deseq2()` runs the full DESeq2 pipeline (DESeqDataSet → DESeq → results)
- **Sensible defaults** — match R package defaults, but allow override via `**kwargs`

## Dependencies

- `rpy2 >= 3.5`
- `pandas >= 1.5`
- `numpy >= 1.23`
- R 4.0+ with BiocManager
