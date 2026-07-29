#include <catch2/catch_test_macros.hpp>
#include <cmath>

#include "golden_reader.hpp"
#include "sctoolbox/hpmap.hpp"

TEST_CASE("hpmap matches MATLAB goldens") {
    auto groups = golden::loadGoldens(std::string(GOLDENS_TEXT_DIR) + "/hplmap_private.gold");
    auto& cases = groups.at("cases_hpmap");
    REQUIRE(!cases.empty());

    for (const auto& c : cases) {
        const auto& zpin = c.in("zp");
        const auto& win = c.in("w");
        const auto& betain = c.in("beta");
        const auto& zin = c.in("z");
        const std::complex<double> cc = c.in("c")(0, 0);
        const auto& qdat = c.in("qdat");

        const int p = static_cast<int>(zpin.rows());
        const int n = static_cast<int>(zin.rows());

        Eigen::VectorXcd zp(p);
        for (int i = 0; i < p; ++i) zp(i) = zpin(i, 0);

        Eigen::VectorXcd w(n), z(n);
        Eigen::VectorXd beta(n);
        for (int i = 0; i < n; ++i) {
            w(i) = win(i, 0);
            z(i) = zin(i, 0);
            beta(i) = betain(i, 0).real();
        }

        const Eigen::VectorXcd wp = sctoolbox::hpmap(zp, w, beta, z, cc, qdat.real());

        const auto& expected = c.out("wp");
        for (int i = 0; i < p; ++i) {
            CHECK(std::abs(wp(i) - expected(i, 0)) < c.tol);
        }
    }
}
