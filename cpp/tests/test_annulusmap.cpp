#include <catch2/catch_test_macros.hpp>
#include <cmath>

#include "golden_reader.hpp"
#include "sctoolbox/annulusmap.hpp"
#include "sctoolbox/annulusmap_internal.hpp"
#include "sctoolbox/polygon.hpp"

namespace {
sctoolbox::Polygon polygonFrom(const golden::Case& c, const std::string& wName, const std::string& aName) {
    const auto& win = c.in(wName);
    const auto& ain = c.in(aName);
    const int n = static_cast<int>(win.rows());
    Eigen::VectorXcd w(n);
    Eigen::VectorXd a(n);
    for (int i = 0; i < n; ++i) {
        w(i) = win(i, 0);
        a(i) = ain(i, 0).real();
    }
    return sctoolbox::Polygon(w, a);
}
}  // namespace

TEST_CASE("qinit produces the right work-array size") {
    auto groups = golden::loadGoldens(std::string(GOLDENS_TEXT_DIR) + "/annulus_private.gold");
    auto& cases = groups.at("cases_qinit");
    REQUIRE(!cases.empty());
    for (const auto& c : cases) {
        sctoolbox::AnnulusData dataz;
        dataz.M = c.inInt("M");
        dataz.N = c.inInt("N");
        const auto& alfa0 = c.in("ALFA0");
        const auto& alfa1 = c.in("ALFA1");
        dataz.ALFA0.resize(dataz.M);
        dataz.ALFA1.resize(dataz.N);
        for (int i = 0; i < dataz.M; ++i) dataz.ALFA0(i) = alfa0(0, i).real();
        for (int i = 0; i < dataz.N; ++i) dataz.ALFA1(i) = alfa1(0, i).real();
        const int nptq = c.inInt("nptq");

        const Eigen::VectorXd qwork = sctoolbox::qinit(dataz, nptq);
        const auto& expectedSize = c.out("qwork_size");
        const int expectedRows = static_cast<int>(expectedSize(0, 0).real());
        const int expectedCols = static_cast<int>(expectedSize(0, 1).real());
        CHECK(static_cast<int>(qwork.size()) == expectedRows * expectedCols);
    }
}

TEST_CASE("AnnulusMap::eval matches MATLAB goldens") {
    auto groups = golden::loadGoldens(std::string(GOLDENS_TEXT_DIR) + "/annulus_private.gold");
    auto& cases = groups.at("cases_anneval");
    REQUIRE(!cases.empty());
    for (const auto& c : cases) {
        const sctoolbox::Polygon outer = polygonFrom(c, "wOuter", "aOuter");
        const sctoolbox::Polygon inner = polygonFrom(c, "wInner", "aInner");
        const sctoolbox::AnnulusMap m(outer, inner);

        const auto& zpin = c.in("zp");
        const auto& expected = c.out("wp");
        const int p = static_cast<int>(zpin.rows());
        for (int i = 0; i < p; ++i) {
            const std::complex<double> wp = m.eval(zpin(i, 0));
            CHECK(std::abs(wp - expected(i, 0)) < c.tol);
        }
    }
}

TEST_CASE("AnnulusMap::evalinv matches MATLAB goldens") {
    auto groups = golden::loadGoldens(std::string(GOLDENS_TEXT_DIR) + "/annulus_private.gold");
    auto& cases = groups.at("cases_annevalinv");
    REQUIRE(!cases.empty());
    for (const auto& c : cases) {
        const sctoolbox::Polygon outer = polygonFrom(c, "wOuter", "aOuter");
        const sctoolbox::Polygon inner = polygonFrom(c, "wInner", "aInner");
        const sctoolbox::AnnulusMap m(outer, inner);

        const auto& wpin = c.in("wp");
        const auto& expected = c.out("zp");
        const int p = static_cast<int>(wpin.rows());
        for (int i = 0; i < p; ++i) {
            const std::complex<double> zp = m.evalinv(wpin(i, 0));
            CHECK(std::abs(zp - expected(i, 0)) < c.tol);
        }
    }
}
