#include <catch2/catch_test_macros.hpp>

#include "golden_reader.hpp"
#include "sctoolbox/gaussj.hpp"

TEST_CASE("gaussj matches MATLAB goldens") {
    auto groups = golden::loadGoldens(std::string(GOLDENS_TEXT_DIR) + "/gaussj.gold");
    auto& cases = groups.at("cases");
    REQUIRE(!cases.empty());

    for (const auto& c : cases) {
        const int n = c.inInt("n");
        const double alf = c.inReal("alf");
        const double bet = c.inReal("bet");

        Eigen::VectorXd z, w;
        sctoolbox::gaussj(n, alf, bet, z, w);

        const auto& zg = c.out("z");
        const auto& wg = c.out("w");
        REQUIRE(z.size() == zg.rows());
        REQUIRE(w.size() == wg.rows());
        for (int i = 0; i < n; ++i) {
            CHECK(std::abs(z(i) - zg(i, 0).real()) < c.tol);
            CHECK(std::abs(w(i) - wg(i, 0).real()) < c.tol);
        }
    }
}
