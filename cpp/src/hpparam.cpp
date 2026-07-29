#include "sctoolbox/hpparam.hpp"

#include <cmath>
#include <limits>
#include <vector>

#include "sctoolbox/hpquad.hpp"
#include "sctoolbox/nesolve.hpp"
#include "sctoolbox/scqdata.hpp"

namespace sctoolbox {

namespace {

// Converts y (length n-3) into the n-1 finite real prevertices, matching
// hpparam.m's/hppfun.m's shared formula:
// cp = cumprod([1;exp(-y)]); z = [0;cumsum(cp)] - [flipud(cumsum(flipud(cp)));0];
// z = z/z(n-1);
Eigen::VectorXd yToZHpFinite(const Eigen::VectorXd& y, int n) {
    const int m = n - 2;  // length of cp
    Eigen::VectorXd cp(m);
    cp(0) = 1.0;
    for (int i = 0; i < n - 3; ++i) cp(i + 1) = std::exp(-y(i));
    double prod = 1.0;
    for (int i = 0; i < m; ++i) {
        prod *= cp(i);
        cp(i) = prod;
    }

    Eigen::VectorXd zp1(n - 1);
    zp1(0) = 0.0;
    double sum = 0.0;
    for (int i = 0; i < m; ++i) {
        sum += cp(i);
        zp1(i + 1) = sum;
    }

    Eigen::VectorXd zp2(n - 1);
    double sum2 = 0.0;
    for (int i = m - 1; i >= 0; --i) {
        sum2 += cp(i);
        zp2(i) = sum2;
    }
    zp2(m) = 0.0;

    Eigen::VectorXd z = zp1 - zp2;
    z /= z(n - 2);
    return z;
}

// Port of hpparam.m's nested hppfun residual function.
Eigen::VectorXd hppfun(const Eigen::VectorXd& y, int n, const Eigen::VectorXd& betaFinite,
                       const Eigen::VectorXd& nmlen, const std::vector<int>& left, const std::vector<int>& right,
                       const std::vector<bool>& cmplx, const Eigen::MatrixXd& qdat) {
    const Eigen::VectorXd zReal = yToZHpFinite(y, n);
    Eigen::VectorXcd z(n - 1);
    for (int i = 0; i < n - 1; ++i) z(i) = zReal(i);

    const int k = static_cast<int>(left.size());
    Eigen::VectorXcd zleft(k), zright(k), mid(k);
    std::vector<int> singL(k), singR(k);
    for (int i = 0; i < k; ++i) {
        zleft(i) = z(left[i]);
        zright(i) = z(right[i]);
        mid(i) = (zleft(i) + zright(i)) / 2.0;
        singL[i] = left[i] + 1;
        singR[i] = right[i] + 1;
    }
    for (int i = 0; i < k; ++i)
        if (cmplx[i]) mid(i) += std::complex<double>(0.0, (zright(i) - zleft(i)).real() / 2.0);

    const Eigen::VectorXcd ints =
        hpquad(zleft, mid, singL, z, betaFinite, qdat) - hpquad(zright, mid, singR, z, betaFinite, qdat);

    std::vector<int> idxF1, idxF2;
    for (int i = 0; i < k; ++i) (cmplx[i] ? idxF2 : idxF1).push_back(i);

    const int n1 = static_cast<int>(idxF1.size());
    const int n2 = static_cast<int>(idxF2.size());

    const double denom1 = std::abs(ints(idxF1[0]));
    const std::complex<double> denom2 = ints(0);

    Eigen::VectorXd F(n - 3);
    int pos = 0;
    for (int j = 1; j < n1; ++j) F(pos++) = std::abs(ints(idxF1[j])) / denom1;
    for (int j = 0; j < n2; ++j) F(pos++) = (ints(idxF2[j]) / denom2).real();
    for (int j = 0; j < n2; ++j) F(pos++) = (ints(idxF2[j]) / denom2).imag();

    F -= nmlen;
    return F;
}

}  // namespace

HpParamResult hpparam(const Eigen::VectorXcd& w, const Eigen::VectorXd& beta, double tol, int method) {
    const int n = static_cast<int>(w.size());
    const int nqpts = std::max(static_cast<int>(std::ceil(-std::log10(tol))), 4);
    const Eigen::VectorXd betaFinite = beta.head(n - 1);
    const Eigen::MatrixXd qdat = scqdata(betaFinite, nqpts);

    std::vector<bool> atinf(n);
    for (int i = 0; i < n; ++i) atinf[i] = beta(i) <= -1.0;

    Eigen::VectorXcd z(n);

    if (n == 3) {
        z(0) = -1.0;
        z(1) = 1.0;
        z(2) = std::complex<double>(std::numeric_limits<double>::infinity(), 0.0);
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

        const Eigen::VectorXd nmlen = nmlenFull.tail(k - 1);

        // Initial guess: z0 = linspace(-1,1,n-1); y0 = log(diff(z0(1:n-2))./diff(z0(2:n-1))).
        Eigen::VectorXd z0(n - 1);
        for (int i = 0; i < n - 1; ++i) z0(i) = -1.0 + 2.0 * i / (n - 2);
        Eigen::VectorXd y0(n - 3);
        for (int i = 0; i < n - 3; ++i) y0(i) = std::log((z0(i + 1) - z0(i)) / (z0(i + 2) - z0(i + 1)));

        const Fvec fvec = [&](const Eigen::VectorXd& y) {
            return hppfun(y, n, betaFinite, nmlen, left, right, cmplx, qdat);
        };

        Eigen::VectorXd details = Eigen::VectorXd::Zero(16);
        details(0) = 0.0;
        details(1) = static_cast<double>(method);
        details(5) = 100.0 * (n - 3);
        details(7) = tol;
        details(8) = std::min(std::pow(std::numeric_limits<double>::epsilon(), 2.0 / 3.0), tol / 10.0);
        details(11) = nqpts;

        const NesolveResult r = nesolve(fvec, y0, details);
        const Eigen::VectorXd zReal = yToZHpFinite(r.xf, n);
        for (int i = 0; i < n - 1; ++i) z(i) = zReal(i);
        z(n - 1) = std::complex<double>(std::numeric_limits<double>::infinity(), 0.0);
    }

    const std::complex<double> mid = (z(0) + z(1)) / 2.0;
    const Eigen::VectorXcd zFinite = z.head(n - 1);
    const Eigen::VectorXd betaFiniteOut = beta.head(n - 1);
    const std::vector<int> singAt2{2}, singAt1{1};
    const Eigen::VectorXcd g1 =
        hpquad(z.segment(1, 1), Eigen::VectorXcd::Constant(1, mid), singAt2, zFinite, betaFiniteOut, qdat);
    const Eigen::VectorXcd g2 =
        hpquad(z.segment(0, 1), Eigen::VectorXcd::Constant(1, mid), singAt1, zFinite, betaFiniteOut, qdat);
    const std::complex<double> g = g1(0) - g2(0);
    const std::complex<double> c = (w(0) - w(1)) / g;

    return HpParamResult{z, c, qdat};
}

}  // namespace sctoolbox
