"""Tests for the sctoolbox Python bindings.

These exercise the binding layer -- that NumPy arrays flow correctly through
every map type, that forward/inverse are consistent, that the analytic
derivative matches a finite difference, and that accessors return sane values.
The heavy numerical correctness (vs. MATLAB) is covered by the C++ golden
tests; here we mostly assert round trips and self-consistency.
"""

import numpy as np
import pytest

import sctoolbox as sc

# Reference hexagon shared by most map types (matches tests/fixturePolygons.m).
HEX6 = np.array([4, 2j, -2 + 4j, -3, -3 - 1j, 2 - 2j], dtype=complex)


def hex_poly():
    return sc.Polygon(HEX6)


# The non-convex L-shape. It subdivides more than once during triangulation,
# exercising the crtriang edge-ordering path that the quad-only C++ goldens
# never covered (see cpp/src/crtriang.cpp).
CRDISK_VERTS = np.array([1j, -1 + 1j, -1 - 1j, 1 - 1j, 1, 0], dtype=complex)


def crdisk_poly():
    return sc.Polygon(CRDISK_VERTS)


# (name, factory, canonical-domain sample points from the MATLAB fixtures)
INTERIOR_MAPS = [
    ("disk", lambda: sc.DiskMap(hex_poly()),
     np.array([0.3 + 0.4j, -0.6 + 0.2j, 0.1 - 0.5j, 0.8 + 0.1j], dtype=complex)),
    ("hpl", lambda: sc.HplMap(hex_poly()),
     np.array([0.5 + 1j, -2 + 0.5j, 3 + 2j, 0 + 1j], dtype=complex)),
    ("strip", lambda: sc.StripMap(hex_poly(), (1, 4)),
     np.array([0 + 0.5j, 1 + 0.3j, -1 + 0.7j, 2 + 0.1j], dtype=complex)),
    ("rect", lambda: sc.RectMap(hex_poly(), (1, 2, 3, 4)),
     np.array([1.5, 1.4 + 3j, -0.6 + 1j, 0.5 + 2j], dtype=complex)),
    ("crdisk", lambda: sc.CrDiskMap(crdisk_poly()),
     np.array([0.3 + 0.4j, -0.6 + 0.2j, 0.1 - 0.5j, 0.5 + 0.1j], dtype=complex)),
]

IDS = [m[0] for m in INTERIOR_MAPS]


@pytest.mark.parametrize("name,factory,zp", INTERIOR_MAPS, ids=IDS)
def test_eval_shape_and_finiteness(name, factory, zp):
    m = factory()
    w = m.eval(zp)
    assert w.shape == zp.shape
    assert np.all(np.isfinite(w))


@pytest.mark.parametrize("name,factory,zp", INTERIOR_MAPS, ids=IDS)
def test_forward_inverse_roundtrip(name, factory, zp):
    m = factory()
    w = m.eval(zp)
    z_back = m.evalinv(w)
    assert np.max(np.abs(z_back - zp)) < 1e-6


@pytest.mark.parametrize("name,factory,zp", INTERIOR_MAPS, ids=IDS)
def test_derivative_matches_finite_difference(name, factory, zp):
    m = factory()
    h = 1e-6
    fd = (m.eval(zp + h) - m.eval(zp - h)) / (2 * h)
    an = m.evaldiff(zp)
    assert np.max(np.abs(fd - an)) < 1e-4


@pytest.mark.parametrize("name,factory,zp", INTERIOR_MAPS, ids=IDS)
def test_accuracy_is_small(name, factory, zp):
    m = factory()
    assert m.accuracy() < 1e-5


def test_out_of_domain_maps_to_nan():
    m = sc.DiskMap(hex_poly())
    w = m.eval(np.array([2.0 + 0j, 0.5 + 0j], dtype=complex))
    assert np.isnan(w[0])            # |z| > 1 is outside the disk
    assert np.isfinite(w[1])


def test_exterior_map_roundtrip():
    m = sc.ExterMap(hex_poly())
    # Interior of the unit disk (canonical domain), avoiding the origin,
    # which maps to infinity.
    zp = np.array([0.3 + 0.4j, -0.6 + 0.2j, 0.1 - 0.5j, 0.8 + 0.1j], dtype=complex)
    w = m.eval(zp)
    assert np.all(np.isfinite(w))
    z_back = m.evalinv(w)
    assert np.max(np.abs(z_back - zp)) < 1e-6


def test_annulus_roundtrip():
    q = np.sqrt(2.0)
    a = 1.0 + q
    outer = sc.Polygon(np.array([a + a * 1j, -a + a * 1j, -a - a * 1j, a - a * 1j],
                                dtype=complex))
    inner = sc.Polygon(np.array([q, q * 1j, -q, -q * 1j], dtype=complex))
    m = sc.AnnulusMap(outer, inner)
    assert 0.0 < m.u < 1.0
    for z in [0.7j, -0.8 + 0j, 0.6 + 0.3j]:
        w = m.eval(complex(z))
        z_back = m.evalinv(complex(w))
        assert abs(z_back - z) < 1e-8


def test_polygon_properties():
    p = hex_poly()
    assert len(p) == 6
    assert not p.is_inf
    # Interior angles (normalized by pi) of a simple n-gon sum to n - 2.
    assert np.isclose(np.sum(p.angle), len(p) - 2)
    assert np.allclose(p.vertex, HEX6) or np.allclose(np.sort_complex(p.vertex),
                                                      np.sort_complex(HEX6))


def test_polygon_explicit_angles():
    verts = np.array([1j, -1 + 1j, -1 - 1j, 1 - 1j, 1, 0], dtype=complex)
    auto = sc.Polygon(verts)
    explicit = sc.Polygon(verts, auto.angle)
    assert np.allclose(auto.angle, explicit.angle)


def test_unsupported_crdisk_polygon_raises_cleanly():
    # This hexagon requires the narrow-channel edge splitting that is not
    # ported to C++. It must raise a catchable RuntimeError rather than
    # crashing the interpreter (regression for the crtriang OOB fix).
    verts = np.array([4, 2j, -2 + 4j, -3, -3 - 1j, 2 - 2j], dtype=complex)
    with pytest.raises(RuntimeError):
        sc.CrDiskMap(sc.Polygon(verts))
    # Interpreter is still healthy: a supported map still works afterwards.
    m = sc.CrDiskMap(crdisk_poly())
    assert m.accuracy() < 1e-5


def test_prevertices_on_unit_circle():
    m = sc.DiskMap(hex_poly())
    z = np.asarray(m.prevertex)
    assert np.allclose(np.abs(z), 1.0, atol=1e-8)
