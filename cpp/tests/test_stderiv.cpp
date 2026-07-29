#include <catch2/catch_test_macros.hpp>
#include <cmath>

#include "golden_reader.hpp"
#include "sctoolbox/stderiv.hpp"

TEST_CASE("stderiv matches MATLAB goldens") {
    auto groups = golden::loadGoldens(std::string(GOLDENS_TEXT_DIR) + "/stripmap_private.gold");
    auto& cases = groups.at("cases_stderiv");
    REQUIRE(!cases.empty());

    for (const auto& c : cases) {
        const auto& zpin = c.in("zp");
        const auto& zin = c.in("z");
        const auto& betain = c.in("beta");
        const std::complex<double> cc = c.in("c")(0, 0);

        const int npts = static_cast<int>(zpin.rows());
        const int n = static_cast<int>(zin.rows());

        Eigen::VectorXcd zp(npts), z(n);
        Eigen::VectorXd beta(n);
        for (int i = 0; i < npts; ++i) zp(i) = zpin(i, 0);
        for (int i = 0; i < n; ++i) {
            z(i) = zin(i, 0);
            beta(i) = betain(i, 0).real();
        }

        const Eigen::VectorXcd fp = sctoolbox::stderiv(zp, z, beta, cc);

        const auto& fpg = c.out("fp");
        REQUIRE(fp.size() == fpg.rows());
        for (int i = 0; i < fpg.rows(); ++i) CHECK(std::abs(fp(i) - fpg(i, 0)) < c.tol);
    }
}
