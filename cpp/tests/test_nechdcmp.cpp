#include <catch2/catch_test_macros.hpp>

#include "golden_reader.hpp"
#include "sctoolbox/nechdcmp.hpp"

TEST_CASE("nechdcmp matches MATLAB goldens") {
    auto groups = golden::loadGoldens(std::string(GOLDENS_TEXT_DIR) + "/nesolve.gold");
    auto& cases = groups.at("cases_nechdcmp");
    REQUIRE(!cases.empty());

    for (const auto& c : cases) {
        const auto& Hin = c.in("H");
        const int n = static_cast<int>(Hin.rows());
        Eigen::MatrixXd H(n, n);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) H(i, j) = Hin(i, j).real();
        const double maxoffl = c.inReal("maxoffl");

        Eigen::MatrixXd L;
        double mu;
        sctoolbox::nechdcmp(H, maxoffl, L, mu);

        const auto& Lg = c.out("L");
        const double mug = c.out("mu")(0, 0).real();
        CHECK(std::abs(mu - mug) < c.tol);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) CHECK(std::abs(L(i, j) - Lg(i, j).real()) < c.tol);
    }
}
