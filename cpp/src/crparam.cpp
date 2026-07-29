#include "sctoolbox/crparam.hpp"

#include <cmath>
#include <vector>

#include "sctoolbox/craffine.hpp"
#include "sctoolbox/crcdt.hpp"
#include "sctoolbox/crembed.hpp"
#include "sctoolbox/crossrat.hpp"
#include "sctoolbox/crqgraph.hpp"
#include "sctoolbox/crquad.hpp"
#include "sctoolbox/crsplit.hpp"
#include "sctoolbox/crtriang.hpp"
#include "sctoolbox/nesolve.hpp"
#include "sctoolbox/scangle.hpp"
#include "sctoolbox/scqdata.hpp"

namespace sctoolbox {

namespace {

// Port of crparam.m's nested crpfun residual function.
Eigen::VectorXd crpfun(const Eigen::VectorXd& x, int n, const Eigen::VectorXd& beta, const Eigen::VectorXcd& crtarget,
                       const QGraph& Q, const Eigen::MatrixXd& qdat) {
    const Eigen::VectorXd crprever = x.array().exp();
    const int n3 = n - 3;
    Eigen::VectorXcd crimage(n3);

    for (int k = 0; k < n3; ++k) {
        const Eigen::VectorXcd prever = crembed(crprever, Q, k);
        const Eigen::Vector4i idx = Q.qlvert.col(k);
        Eigen::VectorXcd z4(4);
        std::vector<int> sing4(4);
        for (int i = 0; i < 4; ++i) {
            z4(i) = prever(idx(i));
            sing4[i] = idx(i) + 1;
        }
        const Eigen::VectorXcd wv = -crquad(z4, sing4, prever, beta, qdat);
        crimage(k) = ((wv(1) - wv(0)) * (wv(3) - wv(2))) / ((wv(2) - wv(1)) * (wv(0) - wv(3)));
    }

    Eigen::VectorXd f(n3);
    for (int k = 0; k < n3; ++k) f(k) = std::log(std::abs(crimage(k) / crtarget(k)));
    return f;
}

}  // namespace

CrParamResult crparam(const Eigen::VectorXcd& wIn, const Eigen::VectorXd& betaIn, double tol, int method) {
    const int nqpts = std::max(static_cast<int>(std::ceil(-std::log10(tol))), 4);

    const CrSplitResult split = crsplit(wIn);
    const Eigen::VectorXcd w = split.w;
    const int n = static_cast<int>(w.size());
    const Eigen::VectorXd beta = scangle(w);

    CrTriangulation tri = crtriang(w);
    tri = crcdt(w, tri);
    const QGraph Q = crqgraph(w, tri);

    const Eigen::MatrixXd qdat = scqdata(beta, nqpts);
    const Eigen::VectorXcd target = crossrat(w, Q);

    const Eigen::VectorXd z0 = target.array().abs().log();

    const Fvec fvec = [&](const Eigen::VectorXd& x) { return crpfun(x, n, beta, target, Q, qdat); };

    Eigen::VectorXd details = Eigen::VectorXd::Zero(16);
    details(0) = 0.0;
    details(1) = static_cast<double>(method);
    details(5) = 100.0 * (n - 3);
    details(7) = tol;
    details(8) = tol / 10.0;
    details(10) = 12.0;  // max step size
    details(11) = nqpts;

    const NesolveResult r = nesolve(fvec, z0, details, /*identityInitialJacobian=*/true);
    const Eigen::VectorXd cr = r.xf.array().exp();

    const Eigen::MatrixXcd aff = craffine(w, beta, cr, Q, tol);

    return CrParamResult{w, beta, cr, aff, Q, split.orig, qdat};
}

}  // namespace sctoolbox
