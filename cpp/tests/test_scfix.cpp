#include <catch2/catch_test_macros.hpp>
#include <cmath>

#include "golden_reader.hpp"
#include "sctoolbox/scfix.hpp"

namespace {
std::vector<int> toAux(const golden::Matrix& m) {
    std::vector<int> v(m.rows());
    for (int i = 0; i < m.rows(); ++i) v[i] = static_cast<int>(std::lround(m(i, 0).real()));
    return v;
}
}  // namespace

TEST_CASE("scfix matches MATLAB goldens") {
    auto groups = golden::loadGoldens(std::string(GOLDENS_TEXT_DIR) + "/scangle_scfix.gold");
    auto& cases = groups.at("cases_scfix");
    REQUIRE(!cases.empty());

    for (const auto& c : cases) {
        const std::string type = c.inStr("type");
        const auto& win = c.in("w");
        const auto& betain = c.in("beta");
        const int n = static_cast<int>(win.rows());

        Eigen::VectorXcd w(n);
        Eigen::VectorXd beta(n);
        for (int i = 0; i < n; ++i) {
            w(i) = win(i, 0);
            beta(i) = betain(i, 0).real();
        }
        std::vector<int> aux = toAux(c.in("aux"));

        sctoolbox::ScfixResult r = sctoolbox::scfix(type, w, beta, aux);

        const auto& wg = c.out("w");
        const auto& bg = c.out("beta");
        const auto& auxg = toAux(c.out("aux"));
        REQUIRE(r.w.size() == wg.rows());
        REQUIRE(r.beta.size() == bg.rows());
        for (int i = 0; i < wg.rows(); ++i) {
            const auto expected = wg(i, 0);
            if (std::isinf(expected.real()) || std::isinf(expected.imag())) {
                CHECK(r.w(i) == expected);
            } else {
                CHECK(std::abs(r.w(i) - expected) < c.tol);
            }
        }
        for (int i = 0; i < bg.rows(); ++i) {
            CHECK(std::abs(r.beta(i) - bg(i, 0).real()) < c.tol);
        }
        REQUIRE(r.aux.size() == auxg.size());
        for (size_t i = 0; i < auxg.size(); ++i) CHECK(r.aux[i] == auxg[i]);
    }
}
