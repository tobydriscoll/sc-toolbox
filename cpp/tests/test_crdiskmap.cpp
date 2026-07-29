#include <catch2/catch_test_macros.hpp>
#include <cmath>

#include "golden_reader.hpp"
#include "sctoolbox/crdiskmap.hpp"

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
