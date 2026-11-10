"""Test-session compatibility shim.

NumPy 2.4 removed ``np.trapz`` (renamed ``np.trapezoid`` in 2.0). Several v0.3-v0.5 scripts and tests call it, and
some of those scripts are pinned by SHA-256 in released provenance manifests, so they are left byte-identical.
For the test session only, restore the old name as an alias of the new function. Run the historical scripts
directly with NumPy < 2.4, or import this alias first.
"""
import numpy as np

if not hasattr(np, 'trapz'):
    np.trapz = np.trapezoid
