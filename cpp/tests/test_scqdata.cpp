#include <catch2/catch_test_macros.hpp>

#include "golden_reader.hpp"
#include "sctoolbox/scqdata.hpp"

TEST_CASE("scqdata matches MATLAB goldens") {
    auto groups = golden::loadGoldens(std::string(GOLDENS_TEXT_DIR) + "/scqdata.gold");
    auto& cases = groups.at("cases");
    REQUIRE(!cases.empty());

    for (const auto& c : cases) {
        const auto& betaIn = c.in("beta");
        Eigen::VectorXd beta(betaIn.rows());
        for (int i = 0; i < beta.size(); ++i) beta(i) = betaIn(i, 0).real();
        const int nqpts = c.inInt("nqpts");

        Eigen::MatrixXd qdat = sctoolbox::scqdata(beta, nqpts);
        const auto& qdatGold = c.out("qdat");

        REQUIRE(qdat.rows() == qdatGold.rows());
        REQUIRE(qdat.cols() == qdatGold.cols());
        for (int r = 0; r < qdat.rows(); ++r) {
            for (int col = 0; col < qdat.cols(); ++col) {
                CHECK(std::abs(qdat(r, col) - qdatGold(r, col).real()) < c.tol);
            }
        }
    }
}
