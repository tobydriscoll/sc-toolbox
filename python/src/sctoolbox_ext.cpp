// nanobind bindings for the SC Toolbox C++ port.
//
// This exposes the seven Schwarz-Christoffel map classes plus the standalone
// AnnulusMap and the Polygon helper. Eigen vectors/matrices convert to and
// from NumPy arrays automatically (nanobind/eigen/dense.h), so the Python API
// is: pass NumPy complex arrays in, get NumPy complex arrays out.

#include <nanobind/nanobind.h>
#include <nanobind/eigen/dense.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/vector.h>

#include "sctoolbox/annulusmap.hpp"
#include "sctoolbox/crdiskmap.hpp"
#include "sctoolbox/diskmap.hpp"
#include "sctoolbox/dinvmap.hpp"
#include "sctoolbox/extermap.hpp"
#include "sctoolbox/hplmap.hpp"
#include "sctoolbox/polygon.hpp"
#include "sctoolbox/rectmap.hpp"
#include "sctoolbox/stripmap.hpp"

namespace nb = nanobind;
using namespace nb::literals;
using namespace sctoolbox;

namespace {

// Bind the eval/evalinv/evaldiff/accuracy + accessor surface shared by the
// six scmap-derived classes (everything except CrDiskMap, whose evalinv
// signature differs, and AnnulusMap, which is standalone).
template <typename MapT, typename Class>
void bind_common(Class& cls) {
    cls.def("eval", &MapT::eval, "zp"_a,
            "Forward map: canonical-domain points -> polygon-interior points. "
            "Points outside the canonical domain map to NaN.")
        .def(
            "evalinv",
            [](const MapT& m, const Eigen::VectorXcd& wp) { return m.evalinv(wp).zp; },
            "wp"_a, "Inverse map: polygon points -> canonical-domain preimages.")
        .def("evaldiff", &MapT::evaldiff, "zp"_a,
             "Derivative f'(zp) of the forward map.")
        .def("accuracy", &MapT::accuracy,
             "Estimated accuracy of the solved map (memoized).")
        .def_prop_ro("polygon", &MapT::polygon)
        .def_prop_ro("prevertex", &MapT::prevertex)
        .def_prop_ro("constant", &MapT::constant)
        .def_prop_ro("qdata", &MapT::qdata);
}

}  // namespace

NB_MODULE(_sctoolbox, m) {
    m.doc() =
        "Schwarz-Christoffel Toolbox -- Python bindings over the C++ port.\n\n"
        "Conformal maps between canonical domains (disk, half-plane, strip, "
        "rectangle) and polygons, plus the doubly connected annulus map.";

    nb::class_<Polygon>(m, "Polygon",
                        "A polygon: complex vertices (counterclockwise) and "
                        "interior angles normalized by pi.")
        .def(nb::init<Eigen::VectorXcd>(), "vertices"_a,
             "Bounded polygon; interior angles are computed from the geometry.")
        .def(nb::init<Eigen::VectorXcd, Eigen::VectorXd>(), "vertices"_a, "angles"_a,
             "Polygon with explicit interior angles (required for unbounded "
             "polygons). Angles are normalized by pi.")
        .def_prop_ro("vertex", &Polygon::vertex)
        .def_prop_ro("angle", &Polygon::angle)
        .def("__len__", &Polygon::length)
        .def_prop_ro("length", &Polygon::length)
        .def_prop_ro("is_inf", &Polygon::isInf,
                     "True if any vertex is at infinity.");

    auto disk = nb::class_<DiskMap>(m, "DiskMap",
                                    "Map from the unit disk to a polygon interior.");
    disk.def(nb::init<Polygon, double>(), "polygon"_a, "tol"_a = 1e-8);
    bind_common<DiskMap>(disk);

    auto hpl = nb::class_<HplMap>(m, "HplMap",
                                  "Map from the upper half-plane to a polygon interior.");
    hpl.def(nb::init<Polygon, double>(), "polygon"_a, "tol"_a = 1e-8);
    bind_common<HplMap>(hpl);

    auto exter = nb::class_<ExterMap>(m, "ExterMap",
                                      "Map from the disk exterior to a polygon exterior.");
    exter.def(nb::init<Polygon, double>(), "polygon"_a, "tol"_a = 1e-8);
    bind_common<ExterMap>(exter);

    auto strip = nb::class_<StripMap>(m, "StripMap",
                                      "Map from the infinite strip to a polygon interior.");
    strip.def(nb::init<Polygon, std::array<int, 2>, double>(), "polygon"_a, "ends"_a,
              "tol"_a = 1e-8,
              "`ends` are the 1-indexed vertex positions mapping to the strip ends.");
    bind_common<StripMap>(strip);

    auto rect = nb::class_<RectMap>(m, "RectMap",
                                    "Map from a rectangle to a generalized quadrilateral.");
    rect.def(nb::init<Polygon, std::array<int, 4>, double>(), "polygon"_a, "corners"_a,
             "tol"_a = 1e-8,
             "`corners` are the 1-indexed positions of the four quadrilateral "
             "corners (first two describe a long side, CCW order).");
    bind_common<RectMap>(rect);
    rect.def("corners", &RectMap::corners,
             "Indices of the rectangle corners within `prevertex`.")
        .def_prop_ro("strip_length", &RectMap::stripL);

    auto crdisk = nb::class_<CrDiskMap>(
        m, "CrDiskMap",
        "Map from the unit disk to a polygon interior, cross-ratio formulation "
        "(robust for elongated/crowded polygons).");
    crdisk.def(nb::init<Polygon, double>(), "polygon"_a, "tol"_a = 1e-8)
        .def("eval", &CrDiskMap::eval, "zp"_a)
        .def("evalinv", &CrDiskMap::evalinv, "wp"_a)
        .def("evaldiff", &CrDiskMap::evaldiff, "zp"_a)
        .def("accuracy", &CrDiskMap::accuracy)
        .def_prop_ro("polygon", &CrDiskMap::polygon)
        .def_prop_ro("crossratio", &CrDiskMap::crossratio)
        .def_prop_ro("qdata", &CrDiskMap::qdata);

    nb::class_<AnnulusMap>(
        m, "AnnulusMap",
        "Map between a canonical annulus and a doubly connected polygonal "
        "region (outer + inner polygon boundary).")
        .def(nb::init<const Polygon&, const Polygon&>(), "outer"_a, "inner"_a)
        .def("eval", &AnnulusMap::eval, "w"_a,
             "Forward map: annulus point -> doubly connected polygon point.")
        .def("evalinv", &AnnulusMap::evalinv, "z"_a,
             "Inverse map: polygon point -> annulus point.")
        .def_prop_ro("M", &AnnulusMap::M)
        .def_prop_ro("N", &AnnulusMap::N)
        .def_prop_ro("u", &AnnulusMap::u,
                     "Inner radius of the canonical annulus.")
        .def_prop_ro("constant", &AnnulusMap::c);
}
