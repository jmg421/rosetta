"""R session management and bidirectional type conversion."""

import numpy as np
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects import numpy2ri, pandas2ri
from rpy2.robjects.conversion import Converter, localconverter

_converter = Converter("rosetta")
_converter += numpy2ri.converter
_converter += pandas2ri.converter
_converter += ro.default_converter


def to_r_matrix(df: pd.DataFrame):
    """Convert pandas DataFrame to R matrix."""
    from ._errors import RDataError
    if not isinstance(df, pd.DataFrame):
        raise RDataError("Expected pandas DataFrame")
    with localconverter(_converter):
        return ro.r["as.matrix"](ro.conversion.get_conversion().py2rpy(df))


def to_r_dataframe(df: pd.DataFrame):
    """Convert pandas DataFrame to R data.frame."""
    from ._errors import RDataError
    if not isinstance(df, pd.DataFrame):
        raise RDataError("Expected pandas DataFrame")
    with localconverter(_converter):
        return ro.conversion.get_conversion().py2rpy(df)


def to_pandas(r_obj) -> pd.DataFrame:
    """Convert R data.frame/matrix to pandas DataFrame."""
    with localconverter(_converter):
        return ro.conversion.get_conversion().rpy2py(r_obj)


def r_call(func_name: str, *args, **kwargs):
    """Call an R function by name."""
    with localconverter(_converter):
        return ro.r[func_name](*args, **kwargs)
