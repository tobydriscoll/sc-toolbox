#include <catch2/catch_test_macros.hpp>
#include <cmath>

#include "golden_reader.hpp"
#include "sctoolbox/crderiv.hpp"

namespace {
sctoolbox::QGraph readQGraph(const golden::Case& c) {
    const auto& qlvertIn = c.in("Qqlvert");
    const auto& qledgeIn = c.in("Qqledge");
    const auto& adjIn = c.in("Qadjacent");

    sctoolbox::QGraph Q;
    Q.qlvert.resize(qlvertIn.rows(), qlvertIn.cols());
    for (int i = 0; i < qlvertIn.rows(); ++i)
        for (int j = 0; j < qlvertIn.cols(); ++j) Q.qlvert(i, j) = static_cast<int>(qlvertIn(i, j).real()) - 1;

    Q.qledge.resize(qledgeIn.rows(), qledgeIn.cols());
    for (int i = 0; i < qledgeIn.rows(); ++i)
        for (int j = 0; j < qledgeIn.cols(); ++j) Q.qledge(i, j) = static_cast<int>(qledgeIn(i, j).real()) - 1;

    Q.adjacent.resize(adjIn.rows(), adjIn.cols());
    for (int i = 0; i < adjIn.rows(); ++i)
        for (int j = 0; j < adjIn.cols(); ++j) Q.adjacent(i, j) = static_cast<int>(adjIn(i, j).real());

    return Q;
}
}  // namespace

TEST_CASE("crderiv matches MATLAB goldens") {
    auto groups = golden::loadGoldens(std::string(GOLDENS_TEXT_DIR) + "/crdiskmap_private.gold");
    auto& cases = groups.at("cases_crderiv");
    REQUIRE(!cases.empty());

    for (const auto& c : cases) {
        const auto& zpin = c.in("zp");
        const auto& betain = c.in("beta");
        const auto& crin = c.in("cr");
        const auto& affin = c.in("aff");
        const auto& wcfixin = c.in("wcfix");

        const int m = static_cast<int>(zpin.rows());
        const int n3 = static_cast<int>(crin.rows());

        Eigen::VectorXcd zp(m);
        for (int i = 0; i < m; ++i) zp(i) = zpin(i, 0);

        const int n = static_cast<int>(betain.rows());
        Eigen::VectorXd beta(n);
        for (int i = 0; i < n; ++i) beta(i) = betain(i, 0).real();

        Eigen::VectorXd cr(n3);
        for (int i = 0; i < n3; ++i) cr(i) = crin(i, 0).real();

        Eigen::MatrixXcd aff(affin.rows(), affin.cols());
        for (int i = 0; i < affin.rows(); ++i)
            for (int j = 0; j < affin.cols(); ++j) aff(i, j) = affin(i, j);

        Eigen::VectorXcd wcfix(5);
        for (int i = 0; i < 5; ++i) wcfix(i) = wcfixin(0, i);

        const sctoolbox::QGraph Q = readQGraph(c);

        const auto result = sctoolbox::crderiv(zp, beta, cr, aff, wcfix, Q);

        const auto& expected = c.out("fp");
        for (int i = 0; i < m; ++i) {
            CHECK(std::abs(result(i) - expected(i, 0)) < c.tol);
        }
    }
}
