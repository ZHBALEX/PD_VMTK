"""Compatibility entry point for the web centerline workflow.

Run this file to open the Three.js web app:

    python Process_VMTK/Process_VMTK_v3.py

Existing code can still import ``CenterlineExtraction`` from this module.
"""

from pathlib import Path
import sys

if __package__ in {None, ""}:
    sys.path.append(str(Path(__file__).resolve().parents[1]))

from Process_VMTK.server import run


def main() -> None:
    run()


def __getattr__(name):
    if name == "CenterlineExtraction":
        from Process_VMTK.viewer import CenterlineExtraction

        return CenterlineExtraction
    raise AttributeError(name)


if __name__ == "__main__":
    main()
