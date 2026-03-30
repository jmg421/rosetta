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

The core v0.1 infrastructure and basic wrappers already exist (bridge layer, type conversion, 6 working wrappers with 33 tests). The contributor will deepen and harden the library over 350 hours:

### Deepen wrappers with publication-critical parameters (Weeks 1–4)

- **DESeq2**: Add `lfcShrink()` support (apeglm/ashr/normal), contrast specification via `results(contrast=...)`, `resultsNames()`, independent filtering control, LFC thresholds
- **edgeR**: Add contrast matrix support, TREAT threshold testing, quasi-likelihood vs likelihood ratio choice
- **limma-voom**: Add `contrast.matrix`, `treat()`, `decideTests()` support
- **clusterProfiler**: Add GSEA (not just ORA), KEGG/Reactome pathway support, custom gene sets
- **Seurat**: Add `FindMarkers()`, `SCTransform` normalization, integration workflows; refactor to builder/pipeline pattern for multi-step analyses
- **phyloseq**: Add ordination methods, differential abundance, filtering/transformation

### Three-tier API architecture (Weeks 5–6)

- **Tier 1**: Quick defaults (current single-call API for 80% of use cases)
- **Tier 2**: Granular kwargs exposing R parameters directly (shrinkage type, contrasts, filtering)
- **Tier 3**: R escape hatch — expose the underlying rpy2 objects for advanced users

### Subprocess+Rscript fallback (Weeks 7–8)

- Implement alternative backend using `subprocess` + `Rscript` + JSON serialization as fallback when rpy2 installation fails
- Automatic backend detection and switching
- Strict version pinning for rpy2 with defensive error handling

### Testing and documentation (Weeks 9–10)

- Expand pytest suite to cover all new parameters with both mocked and real R integration tests
- Validate wrapper output matches direct R output for published datasets (e.g. airway for DESeq2)
- Docstrings with usage examples for all public functions
- Vignette-style tutorial notebook comparing rosetta output to equivalent R code

### Functions (expanded)

| Function | R Package | Key Parameters Added |
|----------|-----------|---------------------|
| `rb.deseq2()` | DESeq2 | contrast, shrinkage, lfc_threshold, alpha |
| `rb.edger()` | edgeR | contrast, test_type, treat_lfc |
| `rb.limma_voom()` | limma | contrast, treat_lfc, decide_tests |
| `rb.enrichment()` | clusterProfiler | method (ORA/GSEA), database (GO/KEGG/Reactome) |
| `rb.seurat()` | Seurat | Builder pattern, FindMarkers, SCTransform |
| `rb.phyloseq()` | phyloseq | ordination, differential abundance |

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
- CO-MENTOR: Matias Salibian-Barrera <matias@stat.ubc.ca> is a Professor of Statistics at the University of British Columbia. Expert R programmer, co-author of the RobStatTM package, and prior GSoC co-supervisor.

## Tests

- **Easy**: Install DESeq2 in R and run the standard vignette example (airway dataset). Report the top 5 differentially expressed genes by adjusted p-value.
- **Medium**: Using rpy2 in Python, replicate the same DESeq2 analysis from Python. Convert the R results to a pandas DataFrame. Post your script and output.
- **Hard**: Write a Python function `deseq2(counts, metadata, design)` that takes pandas DataFrames, runs DESeq2 via rpy2, and returns a pandas DataFrame of results. Include input validation (non-negative integer counts, matching sample names) and at least 3 pytest tests. Post as a GitHub repo or gist.

## Solutions of tests

Contributors, please post a link to your test results here.
