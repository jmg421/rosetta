"""rosetta — Seamless Python wrappers for R bioinformatics packages."""

from ._errors import RDataError, RFormulaError, RPackageMissing
from .wrappers.deseq2 import deseq2
from .wrappers.edger import edger

__all__ = ["deseq2", "edger", "RDataError", "RFormulaError", "RPackageMissing"]
