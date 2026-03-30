"""Tests for rosetta._errors."""

from rosetta._errors import RDataError, RFormulaError, RPackageMissing, RosettaError


def test_hierarchy():
    assert issubclass(RPackageMissing, RosettaError)
    assert issubclass(RFormulaError, RosettaError)
    assert issubclass(RDataError, RosettaError)


def test_rpackagemissing_message():
    exc = RPackageMissing("DESeq2")
    assert "DESeq2" in str(exc)
    assert exc.package == "DESeq2"


def test_rformulaerror():
    exc = RFormulaError("bad formula")
    assert "bad formula" in str(exc)


def test_rdataerror():
    exc = RDataError("negative counts")
    assert "negative counts" in str(exc)
