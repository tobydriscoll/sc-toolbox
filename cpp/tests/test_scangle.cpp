#include <catch2/catch_test_macros.hpp>

#include "golden_reader.hpp"
#include "sctoolbox/scangle.hpp"

TEST_CASE("scangle matches MATLAB goldens") {
    auto groups = golden::loadGoldens(std::string(GOLDENS_TEXT_DIR) + "/scangle_scfix.gold");
    auto& cases = groups.at("cases_angle");
    REQUIRE(!cases.empty());

    for (const auto& c : cases) {
        const auto& w = c.in("w");
        Eigen::VectorXcd wv(w.rows());
        for (int i = 0; i < w.rows(); ++i) wv(i) = w(i, 0);

        Eigen::VectorXd beta = sctoolbox::scangle(wv);
        const auto& betaGold = c.out("beta");

        REQUIRE(beta.size() == betaGold.rows());
        for (int i = 0; i < beta.size(); ++i) {
            const bool nanExpected = std::isnan(betaGold(i, 0).real());
            if (nanExpected) {
                CHECK(std::isnan(beta(i)));
            } else {
                CHECK(std::abs(beta(i) - betaGold(i, 0).real()) < c.tol);
            }
        }
    }
}
