#include <catch2/catch_test_macros.hpp>
#include <cmath>

#include "golden_reader.hpp"
#include "sctoolbox/rderiv.hpp"
#include "sctoolbox/rectmap.hpp"

namespace {
sctoolbox::RectMap buildMap(const golden::Case& c) {
    const auto& win = c.in("w");
    const auto& betain = c.in("beta");
    const int n = static_cast<int>(win.rows());
    Eigen::VectorXcd w(n);
    Eigen::VectorXd beta(n);
    for (int i = 0; i < n; ++i) {
        w(i) = win(i, 0);
        beta(i) = betain(i, 0).real();
    }
    const auto& cornersin = c.in("corners");
    std::array<int, 4> cnr;
    for (int i = 0; i < 4; ++i) cnr[i] = static_cast<int>(cornersin(0, i).real());
    return sctoolbox::RectMap(sctoolbox::Polygon(w, beta.array() + 1.0), cnr, 1e-9);
}
}  // namespace

TEST_CASE("RectMap::eval matches MATLAB goldens") {
    auto groups = golden::loadGoldens(std::string(GOLDENS_TEXT_DIR) + "/rectmap_class.gold");
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

TEST_CASE("RectMap::evalinv matches MATLAB goldens") {
    auto groups = golden::loadGoldens(std::string(GOLDENS_TEXT_DIR) + "/rectmap_class.gold");
    auto& cases = groups.at("cases_evalinv");
    REQUIRE(!cases.empty());
    for (const auto& c : cases) {
        const auto m = buildMap(c);
        const auto& wpin = c.in("wp");
        const int p = static_cast<int>(wpin.rows());
        Eigen::VectorXcd wp(p);
        for (int i = 0; i < p; ++i) wp(i) = wpin(i, 0);
        const auto result = m.evalinv(wp);
        const auto& expected = c.out("zp");
        for (int i = 0; i < p; ++i) CHECK(std::abs(result.zp(i) - expected(i, 0)) < c.tol);
    }
}

TEST_CASE("RectMap::evaldiff matches MATLAB goldens") {
    auto groups = golden::loadGoldens(std::string(GOLDENS_TEXT_DIR) + "/rectmap_class.gold");
    auto& cases = groups.at("cases_evaldiff");
    REQUIRE(!cases.empty());
    for (const auto& c : cases) {
        const auto& zin = c.in("z");
        const auto& betain = c.in("beta");
        const int n = static_cast<int>(zin.rows());
        Eigen::VectorXcd z(n);
        Eigen::VectorXd beta(n);
        for (int i = 0; i < n; ++i) {
            z(i) = zin(i, 0);
            beta(i) = betain(i, 0).real();
        }
        const std::complex<double> cc = c.in("c")(0, 0);
        const double L = c.in("L")(0, 0).real();
        const auto& zpin = c.in("zp");
        const int p = static_cast<int>(zpin.rows());
        Eigen::VectorXcd zp(p);
        for (int i = 0; i < p; ++i) zp(i) = zpin(i, 0);

        const auto fp = sctoolbox::rderiv(zp, z, beta, cc, L);
        const auto& expected = c.out("fp");
        for (int i = 0; i < p; ++i) CHECK(std::abs(fp(i) - expected(i, 0)) < c.tol);
    }
}

TEST_CASE("RectMap::accuracy matches MATLAB goldens") {
    auto groups = golden::loadGoldens(std::string(GOLDENS_TEXT_DIR) + "/rectmap_class.gold");
    auto& cases = groups.at("cases_accuracy");
    REQUIRE(!cases.empty());
    for (const auto& c : cases) {
        const auto m = buildMap(c);
        const double acc = m.accuracy();
        CHECK(std::abs(acc - c.out("acc")(0, 0).real()) < c.tol);
    }
}
