#include <catch2/catch_test_macros.hpp>
#include <cmath>

#include "golden_reader.hpp"
#include "sctoolbox/stderiv.hpp"
#include "sctoolbox/stripmap.hpp"

namespace {
sctoolbox::StripMap buildMap(const golden::Case& c, std::array<int, 2> ends) {
    const auto& win = c.in("w");
    const auto& betain = c.in("beta");
    const int n = static_cast<int>(win.rows());
    Eigen::VectorXcd w(n);
    Eigen::VectorXd beta(n);
    for (int i = 0; i < n; ++i) {
        w(i) = win(i, 0);
        beta(i) = betain(i, 0).real();
    }
    return sctoolbox::StripMap(sctoolbox::Polygon(w, beta.array() + 1.0), ends, 1e-12);
}
}  // namespace

TEST_CASE("StripMap::eval matches MATLAB goldens") {
    auto groups = golden::loadGoldens(std::string(GOLDENS_TEXT_DIR) + "/stripmap_class.gold");
    auto& cases = groups.at("cases_eval");
    REQUIRE(!cases.empty());
    const std::array<int, 2> ends{1, 4};
    for (const auto& c : cases) {
        const auto m = buildMap(c, ends);
        const auto& zpin = c.in("zp");
        const int p = static_cast<int>(zpin.rows());
        Eigen::VectorXcd zp(p);
        for (int i = 0; i < p; ++i) zp(i) = zpin(i, 0);
        const auto wp = m.eval(zp);
        const auto& expected = c.out("wp");
        for (int i = 0; i < p; ++i) CHECK(std::abs(wp(i) - expected(i, 0)) < c.tol);
    }
}

TEST_CASE("StripMap::evalinv matches MATLAB goldens") {
    auto groups = golden::loadGoldens(std::string(GOLDENS_TEXT_DIR) + "/stripmap_class.gold");
    auto& cases = groups.at("cases_evalinv");
    REQUIRE(!cases.empty());
    const std::array<int, 2> ends{1, 4};
    for (const auto& c : cases) {
        const auto m = buildMap(c, ends);
        const auto& wpin = c.in("wp");
        const int p = static_cast<int>(wpin.rows());
        Eigen::VectorXcd wp(p);
        for (int i = 0; i < p; ++i) wp(i) = wpin(i, 0);
        const auto result = m.evalinv(wp);
        const auto& expected = c.out("zp");
        for (int i = 0; i < p; ++i) CHECK(std::abs(result.zp(i) - expected(i, 0)) < c.tol);
    }
}

TEST_CASE("StripMap::evaldiff matches MATLAB goldens") {
    auto groups = golden::loadGoldens(std::string(GOLDENS_TEXT_DIR) + "/stripmap_class.gold");
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
        const auto& zpin = c.in("zp");
        const int p = static_cast<int>(zpin.rows());
        Eigen::VectorXcd zp(p);
        for (int i = 0; i < p; ++i) zp(i) = zpin(i, 0);

        const auto fp = sctoolbox::stderiv(zp, z, beta, cc);
        const auto& expected = c.out("fp");
        for (int i = 0; i < p; ++i) CHECK(std::abs(fp(i) - expected(i, 0)) < c.tol);
    }
}

TEST_CASE("StripMap::accuracy matches MATLAB goldens") {
    auto groups = golden::loadGoldens(std::string(GOLDENS_TEXT_DIR) + "/stripmap_class.gold");
    auto& cases = groups.at("cases_accuracy");
    REQUIRE(!cases.empty());
    const std::array<int, 2> ends{1, 4};
    for (const auto& c : cases) {
        const auto m = buildMap(c, ends);
        const double acc = m.accuracy();
        CHECK(std::abs(acc - c.out("acc")(0, 0).real()) < c.tol);
    }
}
