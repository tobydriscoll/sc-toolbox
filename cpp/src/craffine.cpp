#include "sctoolbox/craffine.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "sctoolbox/crembed.hpp"
#include "sctoolbox/crquad.hpp"
#include "sctoolbox/scqdata.hpp"

namespace sctoolbox {

Eigen::MatrixXcd craffine(const Eigen::VectorXcd& w, const Eigen::VectorXd& beta, const Eigen::VectorXd& cr,
                          const QGraph& Q, double tol) {
    const int n = static_cast<int>(beta.size());
    const int n3 = n - 3;
    const int nqpts = std::max(4, static_cast<int>(std::ceil(-std::log10(tol))));
    const Eigen::MatrixXd qdat = scqdata(beta, nqpts);

    const std::complex<double> nanc(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
    Eigen::VectorXcd wn = Eigen::VectorXcd::Constant(n, nanc);
    Eigen::MatrixXcd aff = Eigen::MatrixXcd::Constant(n3, 2, nanc);
    Eigen::MatrixXcd rawimage = Eigen::MatrixXcd::Ones(4, n3);

    // w is always fully known for crparam's usage, so the first side (0,1)
    // is always usable.
    const int s1 = 0, s2 = 1 % n;
    const std::complex<double> side1 = w(s1), side2 = w(s2);

    int edgenum = -1;
    for (int c = 0; c < Q.edge.cols(); ++c) {
        if (Q.edge(0, c) == s1 && Q.edge(1, c) == s2) {
            edgenum = c;
            break;
        }
    }
    if (edgenum < 0) {
        for (int c = 0; c < Q.edge.cols(); ++c) {
            if (Q.edge(1, c) == s1 && Q.edge(0, c) == s2) {
                edgenum = c;
                break;
            }
        }
    }

    int quadnum = -1;
    for (int q = 0; q < n3; ++q) {
        bool found = false;
        for (int r = 0; r < 4; ++r)
            if (Q.qledge(r, q) == edgenum) found = true;
        if (found) {
            quadnum = q;
            break;
        }
    }

    auto embedAndQuad = [&](int q) {
        const Eigen::VectorXcd z = crembed(cr, Q, q);
        const Eigen::Vector4i idx = Q.qlvert.col(q);
        Eigen::VectorXcd z4(4);
        std::vector<int> sing4(4);
        for (int i = 0; i < 4; ++i) {
            z4(i) = z(idx(i));
            sing4[i] = idx(i) + 1;
        }
        const Eigen::VectorXcd v = -crquad(z4, sing4, z, beta, qdat);
        return std::make_pair(idx, v);
    };

    auto [idx0, v0] = embedAndQuad(quadnum);
    for (int i = 0; i < 4; ++i) wn(idx0(i)) = v0(i);

    Eigen::Matrix2cd A;
    A(0, 0) = wn(s1);
    A(0, 1) = 1.0;
    A(1, 0) = wn(s2);
    A(1, 1) = 1.0;
    Eigen::Vector2cd rhs(side1, side2);
    const Eigen::Vector2cd y0 = A.colPivHouseholderQr().solve(rhs);
    aff(quadnum, 0) = y0(0);
    aff(quadnum, 1) = y0(1);
    for (int i = 0; i < 4; ++i) wn(idx0(i)) = wn(idx0(i)) * aff(quadnum, 0) + aff(quadnum, 1);
    rawimage.col(quadnum) = v0;

    std::vector<int> stack, origin;
    for (int q = 0; q < n3; ++q)
        if (Q.adjacent(q, quadnum) != 0 && std::isnan(aff(q, 0).real())) {
            stack.push_back(q);
            origin.push_back(quadnum);
        }

    while (!stack.empty()) {
        const int q = stack.back();
        stack.pop_back();
        const int oldq = origin.back();
        origin.pop_back();

        auto [idxq, newv] = embedAndQuad(q);
        rawimage.col(q) = newv;

        const Eigen::VectorXcd oldv = rawimage.col(oldq);
        const Eigen::Vector4i oldidx = Q.qlvert.col(oldq);

        std::vector<int> ref, oldref;
        for (int jj = 0; jj < 4; ++jj)
            for (int ii = 0; ii < 4; ++ii)
                if (idxq(ii) == oldidx(jj)) {
                    ref.push_back(ii);
                    oldref.push_back(jj);
                }

        std::vector<bool> done4(4, false);
        for (int r : ref) done4[r] = true;

        const int m = static_cast<int>(ref.size());
        Eigen::MatrixXcd Amat(m, 2);
        Eigen::VectorXcd bvec(m);
        for (int k = 0; k < m; ++k) {
            Amat(k, 0) = newv(ref[k]);
            Amat(k, 1) = 1.0;
            bvec(k) = oldv(oldref[k]);
        }
        const Eigen::Vector2cd y = Amat.colPivHouseholderQr().solve(bvec);

        aff(q, 0) = y(0) * aff(oldq, 0);
        aff(q, 1) = y(1) * aff(oldq, 0) + aff(oldq, 1);
        for (int i = 0; i < 4; ++i)
            if (!done4[i]) wn(idxq(i)) = newv(i) * aff(q, 0) + aff(q, 1);

        for (int qq = 0; qq < n3; ++qq)
            if (Q.adjacent(qq, q) != 0 && std::isnan(aff(qq, 0).real())) {
                stack.push_back(qq);
                origin.push_back(q);
            }
    }

    return aff;
}

}  // namespace sctoolbox
