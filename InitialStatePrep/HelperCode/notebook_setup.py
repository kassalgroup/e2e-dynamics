import sys
import importlib
import importlib.util
from pathlib import Path

# import nbimporter
import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg
from scipy.linalg import ldl

# Directory containing ``notebook_setup.py``.
current_dir = Path(__file__).resolve().parent

# Root of the photodissociation repository.
root_dir = current_dir.parent.parent

# Shared project source directory.
src_dir = root_dir / "src"

if str(src_dir) not in sys.path:
    sys.path.insert(0, str(src_dir))
# Import project modules after making the shared source directory available.
import parse_psi4
import io_functions as io
import nuc_wavefxn
import isp_helper_functions
import physical_consts
import grid_parameters
import isp_resources
import isp_report


# Control the public namespace used by ``from notebook_setup import *``.
__all__ = [
    "np",
    "plt",
    "linalg",
    "ldl",
    "Path",
    "current_dir",
    "root_dir",
    "src_dir",
    "parse_psi4",
    "io",
    "nuc_wavefxn",
    "isp_helper_functions",
    "physical_consts",
    "grid_parameters",
    "isp_resources",
    "isp_report",
]


