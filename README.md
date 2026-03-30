# 🪨 rosetta

**Seamless Python wrappers for popular R bioinformatics packages.**

Stop copy-pasting between R and Python. `rosetta` lets you call DESeq2, edgeR, limma, and other essential R bioinformatics tools directly from Python with pandas DataFrames in, pandas DataFrames out.

## Why?

Bioinformaticians live in two worlds. Your upstream pipeline is Python (preprocessing, ML, visualization), but the gold-standard differential expression tools are in R. The current options:

- **Switch to RStudio mid-analysis** — breaks your workflow
- **rpy2** — powerful but verbose, requires manual R-to-Python type conversion
- **Rewrite in Python** — PyDESeq2 exists but coverage is sparse

`rosetta` wraps the real R packages (not reimplementations) with a clean Python API. You get the exact same statistical methods, just callable from Python.

## Quick Start

```bash
pip install rosetta
```

```python
import rosetta as rb

# DESeq2 differential expression — one function call
results = rb.deseq2(
    counts=counts_df,        # pandas DataFrame (genes × samples)
    metadata=metadata_df,    # pandas DataFrame with conditions
    design="~ condition"
)

# Results come back as a pandas DataFrame
results.head()
#          baseMean  log2FoldChange  lfcSE    pvalue      padj
# GeneA     1523.4           2.31   0.41  3.2e-08   1.1e-06
# GeneB      892.1          -1.87   0.38  5.6e-07   8.4e-06
```

## Supported Packages

| R Package | Status | Python Function |
|-----------|--------|----------------|
| DESeq2 | 🚧 In Progress | `rb.deseq2()` |
| edgeR | 📋 Planned | `rb.edger()` |
| limma | 📋 Planned | `rb.limma_voom()` |
| Seurat | 📋 Planned | `rb.seurat()` |
| clusterProfiler | 📋 Planned | `rb.enrichment()` |
| phyloseq | 📋 Planned | `rb.phyloseq()` |

## How It Works

`rosetta` uses `rpy2` under the hood but handles all the ugly parts:

1. **Automatic type conversion** — pandas DataFrames ↔ R data.frames, numpy arrays ↔ R matrices
2. **R dependency management** — installs required R packages automatically via Bioconductor
3. **Clean error messages** — translates cryptic R errors into Python exceptions
4. **No R syntax required** — every parameter is a Python keyword argument

## Requirements

- Python 3.9+
- R 4.0+ (with Bioconductor)
- rpy2

## Contributing

This project is early stage. The highest-impact contributions:

1. **Wrappers for new packages** — pick one from the Planned list and submit a PR
2. **Real-world testing** — try it on your actual data and report what breaks
3. **Documentation** — usage examples from real analyses

See [CONTRIBUTING.md](CONTRIBUTING.md) for details.

## License

MIT

## Acknowledgments

Built on top of [rpy2](https://rpy2.github.io/) and the incredible R/Bioconductor ecosystem. This project wraps — not replaces — the original R packages. All credit for the statistical methods goes to their respective authors.
