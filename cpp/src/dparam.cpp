#include "sctoolbox/dparam.hpp"

#include <cmath>
#include <limits>
#include <vector>

#include "sctoolbox/dabsquad.hpp"
#include "sctoolbox/dquad.hpp"
#include "sctoolbox/nesolve.hpp"
#include "sctoolbox/scqdata.hpp"

namespace sctoolbox {

namespace {

// Converts the unconstrained nesolve variable y (length n-3) into the n
// disk prevertices, matching dparam.m's/dpfun.m's shared z(y) formula:
// cs = cumsum(cumprod([1;exp(-y)])); theta = pi*cs(1:n-3)/cs(end);
// z(1:n-3) = exp(i*theta); z(n-2:n-1) = [-1;-1i]; z(n) = 1 (unset, default).
Eigen::VectorXcd yToZ(const Eigen::VectorXd& y, int n) {
    const int m = n - 2;  // length of cs
    Eigen::VectorXd v(m);
    v(0) = 1.0;
    for (int i = 0; i < n - 3; ++i) v(i + 1) = std::exp(-y(i));

    Eigen::VectorXd cs(m);
    double prod = 1.0;
    for (int i = 0; i < m; ++i) {
        prod *= v(i);
        cs(i) = prod;
    }
    double sum = 0.0;
    for (int i = 0; i < m; ++i) {
        sum += cs(i);
        cs(i) = sum;
    }

    Eigen::VectorXcd z = Eigen::VectorXcd::Ones(n);
    for (int i = 0; i < n - 3; ++i) {
        const double theta = M_PI * cs(i) / cs(m - 1);
        z(i) = std::exp(std::complex<double>(0.0, theta));
    }
    z(n - 3) = -1.0;
    z(n - 2) = std::complex<double>(0.0, -1.0);
    // z(n-1) keeps its default value of 1 from the Ones() initialization,
    // matching dparam.m's `z = ones(n,1)` never being overwritten at z(n).
    return z;
}

// Port of dparam.m's nested dpfun residual function.
Eigen::VectorXd dpfun(const Eigen::VectorXd& y, int n, const Eigen::VectorXd& beta, const Eigen::VectorXd& nmlen,
                      const std::vector<int>& left, const std::vector<int>& right,
                      const std::vector<bool>& cmplxOrig, const Eigen::MatrixXd& qdat) {
    const Eigen::VectorXcd z = yToZ(y, n);
    const int k = static_cast<int>(left.size());

    Eigen::VectorXcd zleft(k), zright(k);
    for (int i = 0; i < k; ++i) {
        zleft(i) = z(left[i]);
        zright(i) = z(right[i]);
    }

    Eigen::VectorXd angl(k);
    for (int i = 0; i < k; ++i) angl(i) = std::arg(zleft(i));

    Eigen::VectorXcd mid(k);
    for (int i = 0; i < k; ++i) {
        const double rat = std::arg(zright(i) / zleft(i));
        const double rem2pi = std::fmod(rat + 2.0 * M_PI, 2.0 * M_PI);
        mid(i) = std::exp(std::complex<double>(0.0, angl(i) + rem2pi / 2.0));
    }

    std::vector<bool> cmplx = cmplxOrig;
    for (int i = 0; i < k; ++i)
        if (cmplx[i]) mid(i) = 0.0;

    bool anyCmplx = false;
    for (bool b : cmplx)
        if (b) anyCmplx = true;
    std::vector<bool> cmplxForInts = cmplx;
    if (anyCmplx) cmplxForInts[0] = true;

    Eigen::VectorXcd ints(k);
    std::vector<int> idxAbs, idxCmplx;
    for (int i = 0; i < k; ++i) (cmplxForInts[i] ? idxCmplx : idxAbs).push_back(i);

    if (!idxAbs.empty()) {
        const int m = static_cast<int>(idxAbs.size());
        Eigen::VectorXcd zl(m), zr(m), midv(m);
        std::vector<int> singL(m), singR(m);
        for (int j = 0; j < m; ++j) {
            const int i = idxAbs[j];
            zl(j) = zleft(i);
            zr(j) = zright(i);
            midv(j) = mid(i);
            singL[j] = left[i] + 1;
            singR[j] = right[i] + 1;
        }
        const Eigen::VectorXd I1 = dabsquad(zl, midv, singL, z, beta, qdat);
        const Eigen::VectorXd I2 = dabsquad(zr, midv, singR, z, beta, qdat);
        for (int j = 0; j < m; ++j) ints(idxAbs[j]) = I1(j) + I2(j);
    }
    if (!idxCmplx.empty()) {
        const int m = static_cast<int>(idxCmplx.size());
        Eigen::VectorXcd zl(m), zr(m), midv(m);
        std::vector<int> singL(m), singR(m);
        for (int j = 0; j < m; ++j) {
            const int i = idxCmplx[j];
            zl(j) = zleft(i);
            zr(j) = zright(i);
            midv(j) = mid(i);
            singL[j] = left[i] + 1;
            singR[j] = right[i] + 1;
        }
        const Eigen::VectorXcd I1 = dquad(zl, midv, singL, z, beta, qdat);
        const Eigen::VectorXcd I2 = dquad(zr, midv, singR, z, beta, qdat);
        for (int j = 0; j < m; ++j) ints(idxCmplx[j]) = I1(j) - I2(j);
    }

    cmplx[0] = false;
    std::vector<int> idxF1, idxF2;
    for (int i = 0; i < k; ++i) (cmplx[i] ? idxF2 : idxF1).push_back(i);

    const std::complex<double> denom = ints(idxF1[0]);
    const int n1 = static_cast<int>(idxF1.size());
    const int n2 = static_cast<int>(idxF2.size());

    Eigen::VectorXd F(n - 3);
    int pos = 0;
    for (int j = 1; j < n1; ++j) F(pos++) = (ints(idxF1[j]) / std::abs(denom)).real();
    for (int j = 0; j < n2; ++j) F(pos++) = (ints(idxF2[j]) / denom).real();
    for (int j = 0; j < n2; ++j) F(pos++) = (ints(idxF2[j]) / denom).imag();

    F -= nmlen;
    return F;
}

}  // namespace

DParamResult dparam(const Eigen::VectorXcd& w, const Eigen::VectorXd& beta, double tol, int method) {
    const int n = static_cast<int>(w.size());
    const int nqpts = std::max(static_cast<int>(std::ceil(-std::log10(tol))), 4);
    const Eigen::MatrixXd qdat = scqdata(beta, nqpts);

    std::vector<bool> atinf(n);
    for (int i = 0; i < n; ++i) atinf[i] = beta(i) <= -1.0;

    Eigen::VectorXcd z(n);

    if (n == 3) {
        z(0) = std::complex<double>(0.0, -1.0);
        z(1) = std::complex<double>(1.0, -1.0) / std::sqrt(2.0);
        z(2) = 1.0;
    } else {
        std::vector<int> left, right;
        for (int i = 0; i <= n - 3; ++i)
            if (!atinf[i]) left.push_back(i);
        for (int i = 1; i <= n - 2; ++i)
            if (!atinf[i]) right.push_back(i);

        const int k = static_cast<int>(left.size());
        std::vector<bool> cmplx(k);
        for (int i = 0; i < k; ++i) cmplx[i] = (right[i] - left[i] == 2);

        Eigen::VectorXcd nmlenC(k);
        for (int i = 0; i < k; ++i) nmlenC(i) = (w(right[i]) - w(left[i])) / (w(1) - w(0));

        std::vector<int> idxAbs, idxCmplx;
        for (int i = 0; i < k; ++i) (cmplx[i] ? idxCmplx : idxAbs).push_back(i);

        Eigen::VectorXd nmlenFull(k);
        int pos = 0;
        for (int i : idxAbs) nmlenFull(pos++) = std::abs(nmlenC(i));
        for (int i : idxCmplx) nmlenFull(pos++) = nmlenC(i).real();
        for (int i : idxCmplx) nmlenFull(pos++) = nmlenC(i).imag();

        Eigen::VectorXd nmlen = nmlenFull.tail(k - 1);  // drop first entry

        const Eigen::VectorXd y0 = Eigen::VectorXd::Zero(n - 3);

        const Fvec fvec = [&](const Eigen::VectorXd& y) { return dpfun(y, n, beta, nmlen, left, right, cmplx, qdat); };

        Eigen::VectorXd details = Eigen::VectorXd::Zero(16);
        details(0) = 0.0;                                                 // trace
        details(1) = static_cast<double>(method);                         // globalization method
        details(5) = 100.0 * (n - 3);                                     // maxiter
        details(7) = tol;                                                 // ftol
        details(8) = std::min(std::pow(std::numeric_limits<double>::epsilon(), 2.0 / 3.0), tol / 10.0);  // steptol
        details(11) = nqpts;

        const NesolveResult r = nesolve(fvec, y0, details);
        z = yToZ(r.xf, n);
    }

    const std::complex<double> mid = (z(0) + z(1)) / 2.0;
    const std::vector<int> singAt1{1}, singAt2{2};
    // dquad(z(2),mid,2,...) - dquad(z(1),mid,1,...) in MATLAB 1-indexing.
    const Eigen::VectorXcd I1 = dquad(z.segment(1, 1), Eigen::VectorXcd::Constant(1, mid), singAt2, z, beta, qdat);
    const Eigen::VectorXcd I0 = dquad(z.segment(0, 1), Eigen::VectorXcd::Constant(1, mid), singAt1, z, beta, qdat);
    const std::complex<double> c = (w(0) - w(1)) / (I1(0) - I0(0));

    return DParamResult{z, c, qdat};
}

}  // namespace sctoolbox
