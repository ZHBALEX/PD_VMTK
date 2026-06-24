"""Tools for VMTK centerline extraction and duct measurements."""

from .config import CenterlineConfig, DisplayConfig

__all__ = ["CenterlineConfig", "DisplayConfig", "CenterlineExtraction"]


def __getattr__(name):
    if name == "CenterlineExtraction":
        from .viewer import CenterlineExtraction

        return CenterlineExtraction
    raise AttributeError(name)
