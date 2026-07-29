#include <catch2/catch_test_macros.hpp>
#include <cmath>

#include "golden_reader.hpp"
#include "sctoolbox/dquad.hpp"

TEST_CASE("dquad matches MATLAB goldens") {
    auto groups = golden::loadGoldens(std::string(GOLDENS_TEXT_DIR) + "/diskmap_private.gold");
    auto& cases = groups.at("cases_dquad");
    REQUIRE(!cases.empty());

    for (const auto& c : cases) {
        const std::complex<double> z1 = c.in("z1")(0, 0);
        const std::complex<double> z2 = c.in("z2")(0, 0);
        const int sing1 = c.inInt("sing1");
        const auto& zin = c.in("z");
        const auto& betain = c.in("beta");
        const int n = static_cast<int>(zin.rows());

        Eigen::VectorXcd z(n);
        Eigen::VectorXd beta(n);
        for (int i = 0; i < n; ++i) {
            z(i) = zin(i, 0);
            beta(i) = betain(i, 0).real();
        }
        const auto& qdat = c.in("qdat");

        Eigen::VectorXcd z1v(1), z2v(1);
        z1v(0) = z1;
        z2v(0) = z2;
        std::vector<int> sing1v{sing1};

        const Eigen::VectorXcd I = sctoolbox::dquad(z1v, z2v, sing1v, z, beta, qdat.real());

        const auto expected = c.out("I")(0, 0);
        CHECK(std::abs(I(0) - expected) < c.tol);
    }
}
