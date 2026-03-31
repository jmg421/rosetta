"""R package detection and installation via BiocManager."""

from rpy2.robjects.conversion import localconverter

from ._bridge import _converter, _get_base
from ._errors import RPackageMissing


def is_installed(package: str) -> bool:
    """Check if an R package is installed."""
    with localconverter(_converter):
        result = _get_base().requireNamespace(package, quietly=True)
        return bool(result[0])


def ensure_installed(package: str) -> None:
    """Ensure an R package is installed, raising RPackageMissing if not."""
    if not is_installed(package):
        raise RPackageMissing(package)
