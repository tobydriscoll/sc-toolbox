"""Schwarz-Christoffel Toolbox -- Python bindings.

Conformal maps between canonical domains (unit disk, upper half-plane,
infinite strip, rectangle) and polygon interiors/exteriors, plus the doubly
connected annulus map. This package wraps the C++ port of the SC Toolbox via
nanobind.

Typical use::

    import numpy as np
    from sctoolbox import Polygon, DiskMap

    # An L-shaped polygon (vertices counterclockwise).
    verts = np.array([0, 2, 2, 1, 1, 0], dtype=complex) + \
            1j * np.array([0, 0, 1, 1, 2, 2], dtype=float)
    m = DiskMap(Polygon(verts))

    # Map disk points into the polygon.
    z = 0.5 * np.exp(1j * np.linspace(0, 2 * np.pi, 100))
    w = m.eval(z)

All ``eval`` / ``evalinv`` / ``evaldiff`` methods take and return NumPy
complex arrays (the annulus map works on scalars).
"""

from ._sctoolbox import (
    AnnulusMap,
    CrDiskMap,
    DiskMap,
    ExterMap,
    HplMap,
    Polygon,
    RectMap,
    StripMap,
)

__all__ = [
    "Polygon",
    "DiskMap",
    "HplMap",
    "ExterMap",
    "StripMap",
    "RectMap",
    "CrDiskMap",
    "AnnulusMap",
]

__version__ = "0.1.0"
