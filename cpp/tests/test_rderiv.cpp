#include <catch2/catch_test_macros.hpp>
#include <cmath>

#include "golden_reader.hpp"
#include "sctoolbox/rderiv.hpp"

TEST_CASE("rderiv matches MATLAB goldens") {
    auto groups = golden::loadGoldens(std::string(GOLDENS_TEXT_DIR) + "/rectmap_private.gold");
    auto& cases = groups.at("cases_rderiv");
    REQUIRE(!cases.empty());

    for (const auto& c : cases) {
        const auto& zpin = c.in("zp");
        const auto& zin = c.in("z");
        const auto& betain = c.in("beta");
        const std::complex<double> cc = c.in("c")(0, 0);
        const double L = c.in("L")(0, 0).real();

        const int p = static_cast<int>(zpin.rows());
        const int n = static_cast<int>(zin.rows());

        Eigen::VectorXcd zp(p);
        for (int i = 0; i < p; ++i) zp(i) = zpin(i, 0);

        Eigen::VectorXcd z(n);
        Eigen::VectorXd beta(n);
        for (int i = 0; i < n; ++i) {
            z(i) = zin(i, 0);
            beta(i) = betain(i, 0).real();
        }

        const auto result = sctoolbox::rderiv(zp, z, beta, cc, L);

        const auto& expected = c.out("fp");
        for (int i = 0; i < p; ++i) {
            CHECK(std::abs(result(i) - expected(i, 0)) < c.tol);
        }
    }
}
