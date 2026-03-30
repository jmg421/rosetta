## Background

Bioinformaticians routinely switch between Python and R mid-analysis. Upstream pipelines (preprocessing, ML, visualization) are Python, but the gold-standard differential expression tools — DESeq2, edgeR, limma — are R packages. Current options are painful: switch to RStudio mid-workflow, write verbose rpy2 boilerplate with manual type conversion, or use incomplete Python reimplementations (e.g. PyDESeq2).

**rosetta** wraps the real R/Bioconductor packages via rpy2 with a clean Python API: pandas DataFrames in, pandas DataFrames out, one function call per analysis.

```python
import rosetta as rb

results = rb.deseq2(
    counts=counts_df,        # pandas DataFrame (genes × samples)
    metadata=metadata_df,    # pandas DataFrame with conditions
    design="~ condition"
)
# Returns pandas DataFrame with baseMean, log2FoldChange, lfcSE, pvalue, padj
```

## Related work

- **rpy2** — powerful but requires manual R↔Python type conversion and verbose R syntax. rosetta builds on rpy2 but hides the complexity.
- **PyDESeq2** — Python reimplementation of DESeq2. Covers only one package, and reimplementing statistical methods risks subtle correctness differences. rosetta wraps the original R code.
- **diffexpr** — thin rpy2 wrapper for DESeq2 only. Unmaintained, no edgeR/limma support, no dependency management.
- **robStatTM Python wrappers (GSoC 2026)** — same rpy2 wrapping approach for robust statistics. rosetta applies this pattern to the Bioconductor ecosystem.

None of these provide a unified, maintained Python interface across multiple Bioconductor differential expression packages with automatic type conversion and R dependency management.

## Details of your coding project

The contributor will build the following over 350 hours:

### Core infrastructure (Weeks 1–3)

- **`rosetta/_bridge.py`** — R session management (lazy init) and bidirectional type conversion: pandas DataFrame ↔ R data.frame, numpy ndarray ↔ R matrix, Python dict ↔ R named list, None ↔ NULL
- **`rosetta/_deps.py`** — R package detection and automatic installation via `BiocManager::install()` with user confirmation
- **`rosetta/_errors.py`** — Exception classes translating R errors to Python: `RPackageMissing`, `RFormulaError`, `RDataError`

### Wrappers (Weeks 4–8)

- **`rosetta.deseq2()`** — Full DESeq2 pipeline (DESeqDataSetFromMatrix → DESeq → results) in one call
- **`rosetta.edger()`** — edgeR quasi-likelihood pipeline (DGEList → calcNormFactors → estimateDisp → glmQLFit → glmQLFTest)
- **`rosetta.limma_voom()`** — limma-voom pipeline (voom → lmFit → eBayes → topTable)

Each wrapper: validates inputs, calls `ensure_installed()`, converts types, runs the R pipeline, returns a pandas DataFrame.

### Testing and documentation (Weeks 9–10)

- pytest suite with mocked rpy2 for CI environments without R, plus integration tests against real R for environments with R installed
- Docstrings with usage examples for all public functions
- Vignette-style tutorial notebook comparing rosetta output to equivalent R code

### Functions

| Function | R Package | Input | Output |
|----------|-----------|-------|--------|
| `rb.deseq2()` | DESeq2 | counts df, metadata df, design formula | results DataFrame (baseMean, log2FC, pvalue, padj) |
| `rb.edger()` | edgeR | counts df, metadata df, design formula | results DataFrame (logFC, logCPM, F, PValue, FDR) |
| `rb.limma_voom()` | limma | counts df, metadata df, design formula | results DataFrame (logFC, AveExpr, t, P.Value, adj.P.Val) |

## Expected impact

Differential expression analysis is one of the most common tasks in genomics. DESeq2 alone has 30,000+ citations. By providing a clean Python interface to these tools, rosetta:

- Eliminates the R/Python context switch for thousands of bioinformaticians
- Preserves statistical correctness by wrapping (not reimplementing) the original R packages
- Lowers the barrier for Python-first researchers to use gold-standard methods
- Creates a reusable pattern for wrapping additional Bioconductor packages (Seurat, clusterProfiler, phyloseq)

The package will be published on PyPI and maintained by Nodes Bio.

## Mentors

Contributors, please contact mentors below after completing at least one of the tests below.

- EVALUATING MENTOR: John Muirhead-Gould <john@nodes.bio> is the founder of Nodes Bio, Inc., which builds AI-powered biological network visualization tools. Experience with rpy2, Cytoscape.js, and bioinformatics pipelines.
- CO-MENTOR: [TBD — seeking co-mentor with R package development and/or prior GSoC experience]

## Tests

- **Easy**: Install DESeq2 in R and run the standard vignette example (airway dataset). Report the top 5 differentially expressed genes by adjusted p-value.
- **Medium**: Using rpy2 in Python, replicate the same DESeq2 analysis from Python. Convert the R results to a pandas DataFrame. Post your script and output.
- **Hard**: Write a Python function `deseq2(counts, metadata, design)` that takes pandas DataFrames, runs DESeq2 via rpy2, and returns a pandas DataFrame of results. Include input validation (non-negative integer counts, matching sample names) and at least 3 pytest tests. Post as a GitHub repo or gist.

## Solutions of tests

Contributors, please post a link to your test results here.
