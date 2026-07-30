#include <catch2/catch_test_macros.hpp>
#include <cmath>

#include "golden_reader.hpp"
#include "sctoolbox/crdiskmap.hpp"
#include "sctoolbox/polygon.hpp"

namespace {
sctoolbox::CrDiskMap buildMap(const golden::Case& c) {
    const auto& win = c.in("w");
    const auto& betain = c.in("beta");
    const int n = static_cast<int>(win.rows());
    Eigen::VectorXcd w(n);
    Eigen::VectorXd beta(n);
    for (int i = 0; i < n; ++i) {
        w(i) = win(i, 0);
        beta(i) = betain(i, 0).real();
    }
    return sctoolbox::CrDiskMap(sctoolbox::Polygon(w, beta.array() + 1.0), 1e-12);
}
}  // namespace

TEST_CASE("CrDiskMap::eval matches MATLAB goldens") {
    auto groups = golden::loadGoldens(std::string(GOLDENS_TEXT_DIR) + "/crdiskmap_class.gold");
    auto& cases = groups.at("cases_eval");
    REQUIRE(!cases.empty());
    for (const auto& c : cases) {
        const auto m = buildMap(c);
        const auto& zpin = c.in("zp");
        const int p = static_cast<int>(zpin.rows());
        Eigen::VectorXcd zp(p);
        for (int i = 0; i < p; ++i) zp(i) = zpin(i, 0);
        const auto wp = m.eval(zp);
        const auto& expected = c.out("wp");
        for (int i = 0; i < p; ++i) CHECK(std::abs(wp(i) - expected(i, 0)) < c.tol);
    }
}

TEST_CASE("CrDiskMap::evalinv matches MATLAB goldens") {
    auto groups = golden::loadGoldens(std::string(GOLDENS_TEXT_DIR) + "/crdiskmap_class.gold");
    auto& cases = groups.at("cases_evalinv");
    REQUIRE(!cases.empty());
    for (const auto& c : cases) {
        const auto m = buildMap(c);
        const auto& wpin = c.in("wp");
        const int p = static_cast<int>(wpin.rows());
        Eigen::VectorXcd wp(p);
        for (int i = 0; i < p; ++i) wp(i) = wpin(i, 0);
        const auto zp = m.evalinv(wp);
        const auto& expected = c.out("zp");
        for (int i = 0; i < p; ++i) CHECK(std::abs(zp(i) - expected(i, 0)) < c.tol);
    }
}

TEST_CASE("CrDiskMap::evaldiff matches MATLAB goldens") {
    // evaldiff's golden doesn't carry w (crderiv doesn't need it), so build
    // the map from the corresponding eval-group case instead (same
    // polygons, same order) and reuse this group's own beta/cr/aff/wcfix.
    auto groups = golden::loadGoldens(std::string(GOLDENS_TEXT_DIR) + "/crdiskmap_class.gold");
    auto& cases = groups.at("cases_evaldiff");
    auto& evalCases = groups.at("cases_eval");
    REQUIRE(!cases.empty());
    int idx = 0;
    for (const auto& c : cases) {
        const auto m = buildMap(evalCases[idx++]);
        const auto& zpin = c.in("zp");
        const int p = static_cast<int>(zpin.rows());
        Eigen::VectorXcd zp(p);
        for (int i = 0; i < p; ++i) zp(i) = zpin(i, 0);
        const auto fp = m.evaldiff(zp);
        const auto& expected = c.out("fp");
        for (int i = 0; i < p; ++i) CHECK(std::abs(fp(i) - expected(i, 0)) < c.tol);
    }
}

TEST_CASE("CrDiskMap::accuracy matches MATLAB goldens") {
    auto groups = golden::loadGoldens(std::string(GOLDENS_TEXT_DIR) + "/crdiskmap_class.gold");
    auto& cases = groups.at("cases_accuracy");
    REQUIRE(!cases.empty());
    for (const auto& c : cases) {
        const auto m = buildMap(c);
        const double acc = m.accuracy();
        CHECK(std::abs(acc - c.out("acc")(0, 0).real()) < c.tol);
    }
}

// Regression test for the crtriang edge-ordering bug (cpp/src/crtriang.cpp):
// polygons that subdivide more than once (n >= 5) used to store a diagonal
// edge with unsorted endpoints, the base-case edge lookup then missed it, and
// the resulting -1 edge index caused an out-of-bounds Eigen write. The MATLAB
// golden cases only cover quadrilaterals (a single subdivision), so this path
// was never exercised. No MATLAB golden is needed here: we only assert
// self-consistency (forward/inverse round trip) and a small accuracy estimate.
namespace {
double crRoundtripError(const sctoolbox::Polygon& poly) {
    const sctoolbox::CrDiskMap m(poly, 1e-10);
    Eigen::VectorXcd zp(24);
    for (int i = 0; i < 24; ++i)
        zp(i) = 0.5 * std::exp(std::complex<double>(0.0, 2.0 * M_PI * i / 24.0));
    const Eigen::VectorXcd wp = m.eval(zp);
    const Eigen::VectorXcd back = m.evalinv(wp);
    double err = 0.0;
    for (int i = 0; i < 24; ++i) err = std::max(err, std::abs(back(i) - zp(i)));
    return err;
}
}  // namespace

TEST_CASE("CrDiskMap handles polygons that subdivide more than once (n>=5)") {
    // A regular octagon: several subdivisions, exercises the wrap-around
    // sub-polygon path that mis-ordered the stored diagonal edges.
    Eigen::VectorXcd oct(8);
    for (int i = 0; i < 8; ++i)
        oct(i) = std::exp(std::complex<double>(0.0, 2.0 * M_PI * i / 8.0));
    CHECK(crRoundtripError(sctoolbox::Polygon(oct)) < 1e-6);

    // The L-shape from tests/fixturePolygons.m -- a non-convex hexagon.
    Eigen::VectorXcd L(6);
    L << std::complex<double>(0, 1), std::complex<double>(-1, 1),
        std::complex<double>(-1, -1), std::complex<double>(1, -1),
        std::complex<double>(1, 0), std::complex<double>(0, 0);
    CHECK(crRoundtripError(sctoolbox::Polygon(L)) < 1e-6);
}
