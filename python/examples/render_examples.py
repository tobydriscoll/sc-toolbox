"""Render Schwarz-Christoffel maps to PNGs for the README.

Each figure is a two-panel "carpet plot": the left panel shows an orthogonal
grid drawn in the map's canonical domain (disk, half-plane, strip, rectangle,
or annulus); the right panel shows the conformal image of that same grid inside
the target polygon. Because the maps are conformal, the two families of grid
lines stay orthogonal after mapping.

It also renders `gallery.png`, a single six-tile montage (one map type per
tile, grid lines colored by a cyclic colormap) used as the hero image in the
top-level README.

Run from anywhere:

    python examples/render_examples.py            # -> examples/images/*.png

Requires: numpy, matplotlib, and the built `sctoolbox` extension on the path.
"""

from __future__ import annotations

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import sctoolbox as sc

HERE = os.path.dirname(os.path.abspath(__file__))
IMG_DIR = os.path.join(HERE, "images")
os.makedirs(IMG_DIR, exist_ok=True)

# Two grid-line families get two colors; they remain orthogonal under the map.
COLOR_A = "#1f77b4"  # blue  (radial / vertical / spokes)
COLOR_B = "#d62728"  # red   (circular / horizontal)
BOUNDARY = "#222222"

# The reference hexagon used across most interior/exterior examples.
HEX6 = np.array([4, 2j, -2 + 4j, -3, -3 - 1j, 2 - 2j], dtype=complex)


def _finite(z):
    return z[np.isfinite(z)]


def _draw_polygon(ax, poly, **kw):
    """Draw a (bounded) polygon boundary, closing the loop."""
    v = np.asarray(poly.vertex)
    v = _finite(v)
    v = np.append(v, v[0])
    ax.plot(v.real, v.imag, color=BOUNDARY, lw=1.8, **kw)


def _plot_lines(ax, lines, color):
    for ln in lines:
        ax.plot(ln.real, ln.imag, color=color, lw=0.7)


def _finish(ax, title):
    ax.set_aspect("equal", "box")
    ax.set_title(title, fontsize=11)
    ax.tick_params(labelsize=7)


# --------------------------------------------------------------------------
# Canonical-domain grid generators. Each returns (family_a, family_b): two
# lists of complex polylines.
# --------------------------------------------------------------------------
def disk_grid(rmin=0.0, rmax=0.999, n_spoke=24, n_circle=8, npts=400):
    theta = np.linspace(0, 2 * np.pi, npts)
    radii = np.linspace(max(rmin, 1e-3), rmax, n_circle + 1)
    spokes = [np.linspace(rmin, rmax, npts) * np.exp(1j * a)
              for a in np.linspace(0, 2 * np.pi, n_spoke, endpoint=False)]
    circles = [r * np.exp(1j * theta) for r in radii]
    return spokes, circles


def annulus_grid(u, n_spoke=25, n_circle=6, npts=400):
    # Offset the spoke angles so none aligns exactly with a polygon vertex
    # direction (0/90/180/270 deg for these concentric squares), and keep the
    # circles strictly interior -- rays/points through a prevertex hit an
    # unported quadrature branch. The boundaries are drawn separately.
    theta = np.linspace(0, 2 * np.pi, npts)
    radii = np.linspace(u, 1.0, n_circle + 2)[1:-1]
    offset = 0.131
    spokes = [np.linspace(u, 0.999, npts) * np.exp(1j * (offset + a))
              for a in np.linspace(0, 2 * np.pi, n_spoke, endpoint=False)]
    circles = [r * np.exp(1j * theta) for r in radii]
    return spokes, circles


def box_grid(xlim, ylim, n_vert=25, n_horiz=13, npts=400):
    xs = np.linspace(xlim[0], xlim[1], n_vert)
    ys = np.linspace(ylim[0], ylim[1], n_horiz)
    vert = [xs_i + 1j * np.linspace(ylim[0], ylim[1], npts) for xs_i in xs]
    horiz = [np.linspace(xlim[0], xlim[1], npts) + 1j * ys_i for ys_i in ys]
    return vert, horiz


def render_two_panel(fname, title, src_a, src_b, img_a, img_b, poly,
                     src_boundary=None, src_lim=None, img_lim=None):
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(9.5, 4.6))

    _plot_lines(axL, src_a, COLOR_A)
    _plot_lines(axL, src_b, COLOR_B)
    if src_boundary is not None:
        axL.plot(src_boundary.real, src_boundary.imag, color=BOUNDARY, lw=1.8)
    _finish(axL, "canonical domain")
    if src_lim:
        axL.set_xlim(src_lim[0])
        axL.set_ylim(src_lim[1])

    _plot_lines(axR, img_a, COLOR_A)
    _plot_lines(axR, img_b, COLOR_B)
    _draw_polygon(axR, poly)
    _finish(axR, "polygon (image)")
    if img_lim:
        axR.set_xlim(img_lim[0])
        axR.set_ylim(img_lim[1])

    fig.suptitle(title, fontsize=13, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = os.path.join(IMG_DIR, fname)
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("wrote", os.path.relpath(out, HERE))


def mapped(m, lines):
    return [m.eval(ln) for ln in lines]


# --------------------------------------------------------------------------
# Examples
# --------------------------------------------------------------------------
def example_disk():
    m = sc.DiskMap(sc.Polygon(HEX6))
    a, b = disk_grid()
    unit = np.exp(1j * np.linspace(0, 2 * np.pi, 400))
    render_two_panel("disk.png", "DiskMap:  unit disk -> polygon interior",
                     a, b, mapped(m, a), mapped(m, b), m.polygon,
                     src_boundary=unit)


def example_halfplane():
    m = sc.HplMap(sc.Polygon(HEX6))
    z = _finite(np.asarray(m.prevertex))
    lo, hi = z.real.min() - 2, z.real.max() + 2
    height = hi - lo
    a, b = box_grid((lo, hi), (0.0, height), n_vert=27, n_horiz=13)
    render_two_panel("halfplane.png",
                     "HplMap:  upper half-plane -> polygon interior",
                     a, b, mapped(m, a), mapped(m, b), m.polygon,
                     src_lim=((lo, hi), (0, height)))


def example_strip():
    m = sc.StripMap(sc.Polygon(HEX6), (1, 4))
    z = _finite(np.asarray(m.prevertex))
    lo, hi = z.real.min() - 1, z.real.max() + 1
    a, b = box_grid((lo, hi), (0.0, 1.0), n_vert=33, n_horiz=9)
    render_two_panel("strip.png",
                     "StripMap:  infinite strip -> polygon interior",
                     a, b, mapped(m, a), mapped(m, b), m.polygon,
                     src_lim=((lo, hi), (-0.2, 1.2)))


def example_rect():
    m = sc.RectMap(sc.Polygon(HEX6), (1, 2, 3, 4))
    z = _finite(np.asarray(m.prevertex))
    xlim = (z.real.min(), z.real.max())
    ylim = (z.imag.min(), z.imag.max())
    a, b = box_grid(xlim, ylim, n_vert=21, n_horiz=13)
    corners = np.append(np.asarray(m.prevertex)[list(m.corners())],
                        np.asarray(m.prevertex)[m.corners()[0]])
    render_two_panel("rect.png",
                     "RectMap:  rectangle -> generalized quadrilateral",
                     a, b, mapped(m, a), mapped(m, b), m.polygon,
                     src_boundary=corners)


def example_crdisk():
    # The non-convex L-shape, solved with the cross-ratio formulation, which
    # stays well-conditioned where the plain disk map would crowd prevertices.
    verts = np.array([1j, -1 + 1j, -1 - 1j, 1 - 1j, 1, 0], dtype=complex)
    m = sc.CrDiskMap(sc.Polygon(verts))
    a, b = disk_grid(rmax=0.995)
    unit = np.exp(1j * np.linspace(0, 2 * np.pi, 400))
    render_two_panel("crdisk.png",
                     "CrDiskMap:  disk -> non-convex L-shape (cross-ratio)",
                     a, b, mapped(m, a), mapped(m, b), m.polygon,
                     src_boundary=unit)


def example_exterior():
    m = sc.ExterMap(sc.Polygon(HEX6))
    # Source is the unit disk; the center maps to infinity, so avoid r=0.
    a, b = disk_grid(rmin=0.16, rmax=0.999, n_circle=7)
    unit = np.exp(1j * np.linspace(0, 2 * np.pi, 400))
    v = _finite(np.asarray(m.polygon.vertex))
    pad = 4.0
    render_two_panel("exterior.png",
                     "ExterMap:  disk -> polygon exterior",
                     a, b, mapped(m, a), mapped(m, b), m.polygon,
                     src_boundary=unit,
                     img_lim=((v.real.min() - pad, v.real.max() + pad),
                              (v.imag.min() - pad, v.imag.max() + pad)))


def example_annulus():
    q = np.sqrt(2.0)
    a_ = 1.0 + q
    outer = sc.Polygon(np.array([a_ + a_ * 1j, -a_ + a_ * 1j,
                                 -a_ - a_ * 1j, a_ - a_ * 1j], dtype=complex))
    inner = sc.Polygon(np.array([q, q * 1j, -q, -q * 1j], dtype=complex))
    m = sc.AnnulusMap(outer, inner)
    u = m.u

    spokes, circles = annulus_grid(u)

    def ev_one(z):
        try:
            return m.eval(complex(z))
        except RuntimeError:  # ray through a prevertex: unported quad branch
            return complex("nan")

    ev = lambda ln: np.array([ev_one(z) for z in ln])
    img_a = [ev(s) for s in spokes]
    img_b = [ev(c) for c in circles]

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(9.5, 4.6))
    _plot_lines(axL, spokes, COLOR_A)
    _plot_lines(axL, circles, COLOR_B)
    th = np.linspace(0, 2 * np.pi, 400)
    axL.plot(np.cos(th), np.sin(th), color=BOUNDARY, lw=1.8)
    axL.plot(u * np.cos(th), u * np.sin(th), color=BOUNDARY, lw=1.8)
    _finish(axL, "canonical annulus")

    _plot_lines(axR, img_a, COLOR_A)
    _plot_lines(axR, img_b, COLOR_B)
    _draw_polygon(axR, outer)
    _draw_polygon(axR, inner)
    _finish(axR, "doubly connected region (image)")

    fig.suptitle("AnnulusMap:  annulus -> doubly connected polygon",
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = os.path.join(IMG_DIR, "annulus.png")
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("wrote", os.path.relpath(out, HERE))


# --------------------------------------------------------------------------
# A single polished gallery figure for the top-level README: one tile per map
# type, each showing the conformal image of a grid inside the target region.
# The "spoke"/"vertical" family is colored by a cyclic colormap so the
# angle-preserving structure of the maps reads at a glance.
# --------------------------------------------------------------------------
def _carpet_tile(ax, img_a, img_b, polys, title, cmap="twilight", lim=None):
    n = len(img_a)
    colors = plt.get_cmap(cmap)(np.linspace(0, 1, n, endpoint=False))
    for ln, c in zip(img_a, colors):
        ax.plot(ln.real, ln.imag, color=c, lw=0.9)
    for ln in img_b:
        ax.plot(ln.real, ln.imag, color="#33333344", lw=0.6)
    for poly in polys:
        v = _finite(np.asarray(poly.vertex))
        v = np.append(v, v[0])
        ax.plot(v.real, v.imag, color="#111111", lw=2.0)
    ax.set_aspect("equal", "box")
    ax.set_title(title, fontsize=12, fontweight="bold", pad=6)
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    if lim:
        ax.set_xlim(lim[0])
        ax.set_ylim(lim[1])


def example_gallery():
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 8.4))

    # Dense per-polyline sampling so the conformal images stay smooth even
    # where the map stretches a grid line over a long, sharply curved arc.
    NP = 2400

    # 1. Disk -> hexagon
    m = sc.DiskMap(sc.Polygon(HEX6))
    a, b = disk_grid(n_spoke=36, n_circle=10, npts=NP)
    _carpet_tile(axes[0, 0], mapped(m, a), mapped(m, b), [m.polygon],
                 "DiskMap: disk → polygon", cmap="twilight")

    # 2. Rectangle -> generalized quadrilateral
    m = sc.RectMap(sc.Polygon(HEX6), (1, 2, 3, 4))
    z = _finite(np.asarray(m.prevertex))
    a, b = box_grid((z.real.min(), z.real.max()), (z.imag.min(), z.imag.max()),
                    n_vert=28, n_horiz=16, npts=NP)
    _carpet_tile(axes[0, 1], mapped(m, a), mapped(m, b), [m.polygon],
                 "RectMap: rectangle → quadrilateral", cmap="viridis")

    # 3. Strip -> hexagon
    m = sc.StripMap(sc.Polygon(HEX6), (1, 4))
    z = _finite(np.asarray(m.prevertex))
    a, b = box_grid((z.real.min() - 1, z.real.max() + 1), (0.0, 1.0),
                    n_vert=44, n_horiz=11, npts=NP)
    _carpet_tile(axes[0, 2], mapped(m, a), mapped(m, b), [m.polygon],
                 "StripMap: strip → polygon", cmap="plasma")

    # 4. Disk exterior -> polygon exterior
    m = sc.ExterMap(sc.Polygon(HEX6))
    a, b = disk_grid(rmin=0.16, rmax=0.999, n_spoke=36, n_circle=9, npts=NP)
    v = _finite(np.asarray(m.polygon.vertex))
    pad = 4.0
    _carpet_tile(axes[1, 0], mapped(m, a), mapped(m, b), [m.polygon],
                 "ExterMap: disk → polygon exterior", cmap="twilight",
                 lim=((v.real.min() - pad, v.real.max() + pad),
                      (v.imag.min() - pad, v.imag.max() + pad)))

    # 5. Cross-ratio disk -> non-convex L-shape
    Lshape = np.array([1j, -1 + 1j, -1 - 1j, 1 - 1j, 1, 0], dtype=complex)
    m = sc.CrDiskMap(sc.Polygon(Lshape))
    a, b = disk_grid(rmax=0.995, n_spoke=36, n_circle=10, npts=NP)
    _carpet_tile(axes[1, 1], mapped(m, a), mapped(m, b), [m.polygon],
                 "CrDiskMap: disk → L-shape", cmap="viridis")

    # 6. Annulus -> doubly connected region (scalar eval)
    q = np.sqrt(2.0)
    a_ = 1.0 + q
    outer = sc.Polygon(np.array([a_ + a_ * 1j, -a_ + a_ * 1j,
                                 -a_ - a_ * 1j, a_ - a_ * 1j], dtype=complex))
    inner = sc.Polygon(np.array([q, q * 1j, -q, -q * 1j], dtype=complex))
    m = sc.AnnulusMap(outer, inner)
    spokes, circles = annulus_grid(m.u, n_spoke=37, n_circle=7, npts=NP)

    def ev(ln):
        out = np.empty(len(ln), dtype=complex)
        for i, z in enumerate(ln):
            try:
                out[i] = m.eval(complex(z))
            except RuntimeError:
                out[i] = complex("nan")
        return out

    _carpet_tile(axes[1, 2], [ev(s) for s in spokes], [ev(c) for c in circles],
                 [outer, inner], "AnnulusMap: annulus → doubly connected",
                 cmap="plasma")

    fig.suptitle("Schwarz–Christoffel conformal maps  ·  sctoolbox (Python)",
                 fontsize=16, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    out = os.path.join(IMG_DIR, "gallery.png")
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print("wrote", os.path.relpath(out, HERE))


def main():
    example_gallery()
    example_disk()
    example_halfplane()
    example_strip()
    example_rect()
    example_crdisk()
    example_exterior()
    example_annulus()
    print("\nAll figures written to", os.path.relpath(IMG_DIR, HERE))


if __name__ == "__main__":
    main()
