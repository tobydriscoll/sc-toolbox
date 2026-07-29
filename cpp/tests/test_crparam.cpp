#include <catch2/catch_test_macros.hpp>
#include <cmath>

#include "golden_reader.hpp"
#include "sctoolbox/crparam.hpp"

TEST_CASE("crparam matches MATLAB goldens") {
    auto groups = golden::loadGoldens(std::string(GOLDENS_TEXT_DIR) + "/crdiskmap_private.gold");
    auto& cases = groups.at("cases_crparam");
    REQUIRE(!cases.empty());

    for (const auto& c : cases) {
        const auto& win = c.in("w");
        const auto& betain = c.in("beta");

        const int n = static_cast<int>(win.rows());
        Eigen::VectorXcd w(n);
        Eigen::VectorXd beta(n);
        for (int i = 0; i < n; ++i) {
            w(i) = win(i, 0);
            beta(i) = betain(i, 0).real();
        }

        const auto result = sctoolbox::crparam(w, beta, 1e-12);

        const auto& expectedCr = c.out("cr");
        for (int i = 0; i < expectedCr.rows(); ++i) {
            CHECK(std::abs(result.cr(i) - expectedCr(i, 0).real()) < c.tol);
        }
    }
}
