"""rosetta — Seamless Python wrappers for R bioinformatics packages."""

from ._errors import RDataError, RFormulaError, RPackageMissing
from .wrappers.deseq2 import deseq2

__all__ = ["deseq2", "RDataError", "RFormulaError", "RPackageMissing"]
