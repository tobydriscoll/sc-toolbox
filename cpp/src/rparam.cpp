#include "sctoolbox/rparam.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "sctoolbox/ellipkkp.hpp"
#include "sctoolbox/nesolve.hpp"
#include "sctoolbox/r2strip.hpp"
#include "sctoolbox/rectmap_internal.hpp"
#include "sctoolbox/scqdata.hpp"
#include "sctoolbox/stquad.hpp"
#include "sctoolbox/stquadh.hpp"

namespace sctoolbox {

namespace {

// Port of @rectmap/private/rptrnsfm.m. cnr is 1-indexed; cnr[0] is always 1
// after rparam's renumbering. y has length n-3.
Eigen::VectorXcd rptrnsfm(const Eigen::VectorXd& y, const std::array<int, 4>& cnr) {
    const int n = static_cast<int>(y.size()) + 3;
    const int c1 = cnr[0], c2 = cnr[1], c3 = cnr[2], c4 = cnr[3];
    Eigen::VectorXcd z = Eigen::VectorXcd::Zero(n);
    auto Z = [&](int p) -> std::complex<double>& { return z(p - 1); };
    auto Y = [&](int p) { return y(p - 1); };

    {
        double cum = 0.0;
        int yidx = c1;
        for (int p = c1 + 1; p <= c2 - 1; ++p) {
            cum += std::exp(Y(yidx));
            Z(p) = cum;
            ++yidx;
        }
    }
    {
        double cum = 0.0;
        int yidx = c4 - 3;
        for (int p = c4 - 1; p >= c3 + 1; --p) {
            cum += std::exp(Y(yidx));
            Z(p) = std::complex<double>(0.0, 1.0) + cum;
            --yidx;
        }
    }

    const double xr1 = Z(c2 - 1).real();
    const double xr2 = Z(c3 + 1).real();
    const double meanxr = (xr1 + xr2) / 2.0;
    const double halfdiff = xr2 / 2.0 - xr1 / 2.0;
    Z(c2) = meanxr + std::sqrt(halfdiff * halfdiff + std::exp(2.0 * Y(c2 - 1)));
    Z(c3) = std::complex<double>(0.0, 1.0) + Z(c2);
    Z(c4) = std::complex<double>(0.0, 1.0);

    const double epsD = std::numeric_limits<double>::epsilon();

    // Trefethen-style short-edge transform, shared by the two short edges.
    auto trefethen = [&](const std::vector<double>& ys, std::complex<double> zc2) {
        const int L = static_cast<int>(ys.size());
        Eigen::VectorXd cp(L + 1);
        cp(0) = 1.0;
        for (int j = 0; j < L; ++j) cp(j + 1) = std::exp(-ys[j]);
        double prod = 1.0;
        for (int j = 0; j <= L; ++j) {
            prod *= cp(j);
            cp(j) = prod;
        }
        Eigen::VectorXd cs(L + 1);
        double sum = 0.0;
        for (int j = 0; j <= L; ++j) {
            sum += cp(j);
            cs(j) = sum;
        }
        Eigen::VectorXd cpRev(L + 1);
        for (int j = 0; j <= L; ++j) cpRev(j) = cp(L - j);
        Eigen::VectorXd csRev(L + 1);
        double sum2 = 0.0;
        for (int j = 0; j <= L; ++j) {
            sum2 += cpRev(j);
            csRev(j) = sum2;
        }
        Eigen::VectorXd part2(L + 1);
        for (int j = 0; j <= L; ++j) part2(j) = csRev(L - j);

        Eigen::VectorXd xfull(L + 2);
        xfull(0) = 0.0 - part2(0);
        for (int j = 0; j < L; ++j) xfull(j + 1) = cs(j) - part2(j + 1);
        xfull(L + 1) = cs(L) - 0.0;

        const double xend = xfull(L + 1);
        Eigen::VectorXcd u(L);
        for (int j = 0; j < L; ++j) {
            const double xmid = xfull(j + 1) / xend;
            if (std::abs(xmid) < epsD) {
                u(j) = -zc2 / epsD;
            } else {
                u(j) = std::log(std::complex<double>(xmid, 0.0)) / M_PI;
            }
        }
        return u;
    };

    {
        std::vector<double> ys;
        for (int p = c2; p <= c3 - 2; ++p) ys.push_back(Y(p));
        const Eigen::VectorXcd u = trefethen(ys, Z(c2));
        int idx = 0;
        for (int p = c2 + 1; p <= c3 - 1; ++p) {
            Z(p) = std::complex<double>(0.0, u(idx).imag()) + Z(c2).real() - u(idx).real();
            ++idx;
        }
    }
    {
        std::vector<double> ys;
        for (int p = c4 - 2; p <= n - 3; ++p) ys.push_back(Y(p));
        for (int p = 1; p <= c1 - 1; ++p) ys.push_back(Y(p));
        const Eigen::VectorXcd u = trefethen(ys, Z(c2));
        int idx = 0;
        for (int p = c4 + 1; p <= n; ++p) {
            Z(p) = u(idx);
            ++idx;
        }
        for (int p = 1; p <= c1 - 1; ++p) {
            Z(p) = u(idx);
            ++idx;
        }
    }

    return z;
}

// Port of rparam.m's nested rpfun residual function.
Eigen::VectorXd rpfun(const Eigen::VectorXd& y, int n, const Eigen::VectorXd& beta, const Eigen::VectorXcd& nmlen,
                      const std::vector<int>& left, const std::vector<int>& right, const std::vector<bool>& cmplx,
                      const Eigen::MatrixXd& qdat, const std::array<int, 4>& cnr) {
    const Eigen::VectorXcd z = rptrnsfm(y, cnr);
    const int m = static_cast<int>(left.size());

    Eigen::VectorXcd zleft(m), zright(m);
    for (int i = 0; i < m; ++i) {
        zleft(i) = z(left[i] - 1);
        zright(i) = z(right[i] - 1);
    }

    const auto [e0, e1] = rectStripEnds(z);
    const int ends1 = e0 + 1, ends2 = e1 + 1;

    const double inf = std::numeric_limits<double>::infinity();
    Eigen::VectorXcd zAug(n + 2);
    Eigen::VectorXd betaAug(n + 2);
    {
        int k = 0;
        for (int p = 1; p <= ends1; ++p) {
            zAug(k) = z(p - 1);
            betaAug(k) = beta(p - 1);
            ++k;
        }
        zAug(k) = std::complex<double>(inf, 0.0);
        betaAug(k) = 0.0;
        ++k;
        for (int p = ends1 + 1; p <= ends2; ++p) {
            zAug(k) = z(p - 1);
            betaAug(k) = beta(p - 1);
            ++k;
        }
        zAug(k) = std::complex<double>(-inf, 0.0);
        betaAug(k) = 0.0;
        ++k;
        for (int p = ends2 + 1; p <= n; ++p) {
            zAug(k) = z(p - 1);
            betaAug(k) = beta(p - 1);
            ++k;
        }
    }
    const Eigen::MatrixXd qdatAug = rectAugQdat(qdat, n, e0, e1);

    std::vector<int> leftAug(m), rightAug(m);
    for (int i = 0; i < m; ++i) {
        leftAug[i] = left[i] + (left[i] > ends1 ? 1 : 0) + (left[i] > ends2 ? 1 : 0);
        rightAug[i] = right[i] + (right[i] > ends1 ? 1 : 0) + (right[i] > ends2 ? 1 : 0);
    }

    Eigen::VectorXcd ints(m);
    std::vector<int> idsOnSide, idsAcross;
    for (int i = 0; i < m; ++i) {
        const bool s2 = (rightAug[i] - leftAug[i] == 1) && (zleft(i).imag() == zright(i).imag());
        (s2 ? idsOnSide : idsAcross).push_back(i);
    }

    if (!idsOnSide.empty()) {
        const int mm = static_cast<int>(idsOnSide.size());
        Eigen::VectorXcd zlv(mm), zrv(mm), midv(mm);
        std::vector<int> singL(mm), singR(mm);
        for (int j = 0; j < mm; ++j) {
            const int i = idsOnSide[j];
            zlv(j) = zleft(i);
            zrv(j) = zright(i);
            midv(j) = (zleft(i) + zright(i)) / 2.0;
            singL[j] = leftAug[i];
            singR[j] = rightAug[i];
        }
        const Eigen::VectorXcd I1 = stquadh(zlv, midv, singL, zAug, betaAug, qdatAug);
        const Eigen::VectorXcd I2 = stquadh(zrv, midv, singR, zAug, betaAug, qdatAug);
        for (int j = 0; j < mm; ++j) ints(idsOnSide[j]) = I1(j) - I2(j);
    }
    if (!idsAcross.empty()) {
        const int mm = static_cast<int>(idsAcross.size());
        Eigen::VectorXcd zlv(mm), zrv(mm), mid1v(mm), mid2v(mm);
        std::vector<int> singL(mm), singR(mm), zeroSing(mm, 0);
        for (int j = 0; j < mm; ++j) {
            const int i = idsAcross[j];
            zlv(j) = zleft(i);
            zrv(j) = zright(i);
            mid1v(j) = std::complex<double>(zleft(i).real(), 0.5);
            mid2v(j) = std::complex<double>(zright(i).real(), 0.5);
            singL[j] = leftAug[i];
            singR[j] = rightAug[i];
        }
        const Eigen::VectorXcd I1 = stquad(zlv, mid1v, singL, zAug, betaAug, qdatAug);
        const Eigen::VectorXcd I2 = stquadh(mid1v, mid2v, zeroSing, zAug, betaAug, qdatAug);
        const Eigen::VectorXcd I3 = stquad(zrv, mid2v, singR, zAug, betaAug, qdatAug);
        for (int j = 0; j < mm; ++j) ints(idsAcross[j]) = I1(j) + I2(j) - I3(j);
    }

    std::vector<int> idxAbs, idxCmplx;
    for (int i = 0; i < m; ++i) (cmplx[i] ? idxCmplx : idxAbs).push_back(i);
    const int nAbs = static_cast<int>(idxAbs.size());

    Eigen::VectorXd absF(nAbs);
    for (int j = 0; j < nAbs; ++j) absF(j) = std::abs(ints(idxAbs[j]));

    std::vector<bool> cmplx2(cmplx.begin() + 1, cmplx.end());
    std::vector<int> nmIdxAbs, nmIdxCmplx;
    for (int i = 0; i < static_cast<int>(cmplx2.size()); ++i) (cmplx2[i] ? nmIdxCmplx : nmIdxAbs).push_back(i);

    Eigen::VectorXd F1(nAbs - 1);
    for (int j = 1; j < nAbs; ++j) F1(j - 1) = std::log((absF(j) / absF(0)) / nmlen(nmIdxAbs[j - 1]).real());

    bool anyCmplx = false;
    for (bool b : cmplx)
        if (b) anyCmplx = true;

    Eigen::VectorXcd F2;
    if (anyCmplx) {
        const int nC = static_cast<int>(idxCmplx.size());
        F2.resize(nC);
        const std::complex<double> denom2 = ints(0);
        for (int j = 0; j < nC; ++j) F2(j) = std::log((ints(idxCmplx[j]) / denom2) / nmlen(nmIdxCmplx[j]));
    }

    Eigen::VectorXd F(F1.size() + (anyCmplx ? 2 * static_cast<int>(F2.size()) : 0));
    int pos = 0;
    for (int j = 0; j < F1.size(); ++j) F(pos++) = F1(j);
    if (anyCmplx) {
        for (int j = 0; j < F2.size(); ++j) F(pos++) = F2(j).real();
        for (int j = 0; j < F2.size(); ++j) F(pos++) = F2(j).imag();
    }
    return F;
}

}  // namespace

RParamResult rparam(const Eigen::VectorXcd& wIn, const Eigen::VectorXd& betaIn, std::array<int, 4> cnrIn,
                    double tol, int method) {
    const int n = static_cast<int>(wIn.size());

    std::vector<int> renum(n);
    const int start0 = cnrIn[0] - 1;
    for (int i = 0; i < n; ++i) renum[i] = (start0 + i) % n;

    Eigen::VectorXcd w(n);
    Eigen::VectorXd beta(n);
    for (int i = 0; i < n; ++i) {
        w(i) = wIn(renum[i]);
        beta(i) = betaIn(renum[i]);
    }

    std::array<int, 4> cnr;
    for (int j = 0; j < 4; ++j) {
        int v = ((cnrIn[j] - cnrIn[0] + 1 + n - 1) % n) + 1;
        cnr[j] = v;
    }
    const int c1 = cnr[0], c2 = cnr[1], c3 = cnr[2], c4 = cnr[3];

    const int nqpts = std::max(static_cast<int>(std::ceil(-std::log10(tol))), 4);
    const Eigen::MatrixXd qdat = scqdata(beta, nqpts);

    std::vector<bool> atinf(n);
    for (int i = 0; i < n; ++i) atinf[i] = beta(i) <= -1.0;

    // Initial guess (z0 = [] branch).
    Eigen::VectorXd dw(n);
    for (int p = 1; p <= n; ++p) dw(p - 1) = std::abs(w(p % n) - w(p - 1));
    {
        double sum = 0.0;
        int cnt = 0;
        for (int i = 0; i < n; ++i)
            if (!std::isinf(dw(i))) {
                sum += dw(i);
                ++cnt;
            }
        const double meanFinite = (cnt > 0) ? sum / cnt : 0.0;
        for (int i = 0; i < n; ++i)
            if (std::isinf(dw(i))) dw(i) = meanFinite;
    }
    auto sumRange = [&](int a, int b) {  // 1-indexed inclusive [a,b], no wraparound
        double s = 0.0;
        for (int p = a; p <= b; ++p) s += dw(p - 1);
        return s;
    };
    const double len = (sumRange(c1, c2 - 1) + sumRange(c3, c4 - 1)) / 2.0;
    const double wid = (sumRange(c2, c3 - 1) + sumRange(c4, n) + sumRange(1, c1 - 1)) / 2.0;
    const double modest = std::min(len / wid, 100.0);

    Eigen::VectorXcd z0v = Eigen::VectorXcd::Zero(n);
    auto Z0 = [&](int p) -> std::complex<double>& { return z0v(p - 1); };
    {
        const int len_ = c2 - c1 + 1;
        for (int j = 0; j < len_; ++j) Z0(c1 + j) = modest * j / (len_ - 1);
    }
    const double dx1 = Z0(c1 + 1).real() - Z0(c1).real();
    for (int p = c1 - 1; p >= 1; --p) Z0(p) = Z0(c1).real() - dx1 * (c1 - p);
    {
        const int cnt = c3 - c2 - 1;
        for (int j = 1; j <= cnt; ++j) Z0(c2 + j) = Z0(c2).real() + dx1 * j;
    }
    {
        const int len_ = c4 - c3 + 1;
        for (int j = 0; j < len_; ++j) Z0(c4 - j) = std::complex<double>(0.0, 1.0) + modest * j / (len_ - 1);
    }
    const double dx2 = (Z0(c4 - 1) - Z0(c4)).real();
    {
        const int cnt = n - c4;
        for (int j = 1; j <= cnt; ++j) Z0(c4 + j) = Z0(c4) - dx2 * j;
    }

    // Convert z0 to unconstrained y0.
    Eigen::VectorXcd dz(n - 1);
    auto Dz = [&](int p) -> std::complex<double>& { return dz(p - 1); };  // 1-indexed
    for (int p = 1; p <= n - 1; ++p) Dz(p) = Z0(p + 1) - Z0(p);
    for (int p = c3; p <= n - 1; ++p) Dz(p) = -Dz(p);

    Eigen::VectorXd y0 = Eigen::VectorXd::Zero(n - 3);
    auto Y0 = [&](int p) -> double& { return y0(p - 1); };
    for (int p = 1; p <= c2 - 2; ++p) Y0(p) = std::log(Dz(p)).real();
    for (int p = c3 - 1; p <= c4 - 3; ++p) Y0(p) = std::log(Dz(p + 2)).real();
    Y0(c2 - 1) = (std::log(Dz(c2 - 1)).real() + std::log(Dz(c3)).real()) / 2.0;

    const double Lguess = Z0(c2).real() - Z0(c1).real();
    {
        const int cnt = c3 - c2 - 1;
        Eigen::VectorXd xv(cnt);
        for (int j = 0; j < cnt; ++j) {
            const std::complex<double> zc = Z0(c2 + 1 + j);
            xv(j) = std::exp(M_PI * (Lguess - std::conj(zc))).real();
        }
        // dx = -diff([1; x; -1])
        Eigen::VectorXd ext(cnt + 2);
        ext(0) = 1.0;
        for (int j = 0; j < cnt; ++j) ext(j + 1) = xv(j);
        ext(cnt + 1) = -1.0;
        Eigen::VectorXd dxFull(cnt + 1);
        for (int j = 0; j <= cnt; ++j) dxFull(j) = -(ext(j + 1) - ext(j));
        for (int j = 0; j < cnt; ++j) Y0(c2 + j) = std::log(dxFull(j) / dxFull(j + 1));
    }
    {
        const int cnt = n - c4;
        Eigen::VectorXd xv(cnt);
        for (int j = 0; j < cnt; ++j) xv(j) = std::exp(M_PI * Z0(c4 + 1 + j)).real();
        Eigen::VectorXd ext(cnt + 2);
        ext(0) = -1.0;
        for (int j = 0; j < cnt; ++j) ext(j + 1) = xv(j);
        ext(cnt + 1) = 1.0;
        Eigen::VectorXd dxFull(cnt + 1);
        for (int j = 0; j <= cnt; ++j) dxFull(j) = ext(j + 1) - ext(j);
        for (int j = 0; j < cnt; ++j) Y0(c4 - 2 + j) = std::log(dxFull(j) / dxFull(j + 1));
    }

    // left/right/cmplx/nmlen.
    std::vector<int> left, right;
    for (int p = 1; p <= n - 2; ++p) left.push_back(p);
    for (int p = 2; p <= n - 1; ++p) right.push_back(p);
    {
        std::vector<int> leftFiltered;
        for (int v : left)
            if (!atinf[v - 1]) leftFiltered.push_back(v);
        left = leftFiltered;
        std::vector<int> rightFiltered;
        for (int v : right)
            if (!atinf[v - 1]) rightFiltered.push_back(v);
        right = rightFiltered;
    }
    if (atinf[n - 2]) right.push_back(n);

    const int m = static_cast<int>(left.size());
    std::vector<bool> cmplx(m);
    for (int i = 0; i < m; ++i) cmplx[i] = (right[i] - left[i] == 2);
    int cmplxCount = 0;
    for (bool b : cmplx)
        if (b) ++cmplxCount;
    if (static_cast<int>(cmplx.size()) + cmplxCount > n - 2) cmplx.back() = false;

    Eigen::VectorXcd nmlenC(m);
    for (int i = 0; i < m; ++i) {
        nmlenC(i) = (w(right[i] - 1) - w(left[i] - 1)) / (w(1) - w(0));
        if (!cmplx[i]) nmlenC(i) = std::abs(nmlenC(i));
    }
    const Eigen::VectorXcd nmlen = nmlenC.tail(m - 1);

    const Fvec fvec = [&](const Eigen::VectorXd& y) { return rpfun(y, n, beta, nmlen, left, right, cmplx, qdat, cnr); };

    Eigen::VectorXd details = Eigen::VectorXd::Zero(16);
    details(0) = 0.0;
    details(1) = static_cast<double>(method);
    details(5) = 100.0 * (n - 3);
    details(7) = tol;
    details(8) = std::min(std::pow(std::numeric_limits<double>::epsilon(), 2.0 / 3.0), tol / 10.0);
    details(11) = nqpts;

    const NesolveResult r = nesolve(fvec, y0, details);
    const Eigen::VectorXcd zStrip = rptrnsfm(r.xf, cnr);

    const auto [e0, e1] = rectStripEnds(zStrip);
    const int ends1 = e0 + 1, ends2 = e1 + 1;
    const double inf = std::numeric_limits<double>::infinity();
    Eigen::VectorXcd zsAug(n + 2);
    Eigen::VectorXd bsAug(n + 2);
    {
        int k = 0;
        for (int p = 1; p <= ends1; ++p) {
            zsAug(k) = zStrip(p - 1);
            bsAug(k) = beta(p - 1);
            ++k;
        }
        zsAug(k) = std::complex<double>(inf, 0.0);
        bsAug(k) = 0.0;
        ++k;
        for (int p = ends1 + 1; p <= ends2; ++p) {
            zsAug(k) = zStrip(p - 1);
            bsAug(k) = beta(p - 1);
            ++k;
        }
        zsAug(k) = std::complex<double>(-inf, 0.0);
        bsAug(k) = 0.0;
        ++k;
        for (int p = ends2 + 1; p <= n; ++p) {
            zsAug(k) = zStrip(p - 1);
            bsAug(k) = beta(p - 1);
            ++k;
        }
    }
    const Eigen::MatrixXd qsAug = rectAugQdat(qdat, n, e0, e1);
    const std::complex<double> midc = (zsAug(0) + zsAug(1)) / 2.0;
    const std::vector<int> singAt2{2}, singAt1{1};
    const Eigen::VectorXcd g1 = stquad(zsAug.segment(1, 1), Eigen::VectorXcd::Constant(1, midc), singAt2, zsAug, bsAug, qsAug);
    const Eigen::VectorXcd g2 = stquad(zsAug.segment(0, 1), Eigen::VectorXcd::Constant(1, midc), singAt1, zsAug, bsAug, qsAug);
    const std::complex<double> g = g1(0) - g2(0);
    const std::complex<double> c = -(w(1) - w(0)) / g;

    // Find prevertices on the rectangle.
    const double Lval = zStrip(c2 - 1).real();
    const EllipKKp kk = ellipkkp(Lval);
    const double K = kk.K, Kp = kk.Kp;
    const std::array<std::complex<double>, 4> rect = {std::complex<double>(K, 0.0), std::complex<double>(K, Kp),
                                                       std::complex<double>(-K, Kp), std::complex<double>(-K, 0.0)};

    std::vector<bool> lSide(n, false), rSide(n, false), tl(n, false), tr(n, false), bl(n, false), brr(n, false);
    for (int p = c3; p <= c4; ++p) lSide[p - 1] = true;
    for (int p = c1; p <= c2; ++p) rSide[p - 1] = true;
    for (int i = 0; i < n; ++i) {
        const double re = zStrip(i).real(), im = zStrip(i).imag();
        if (re > Lval && im > 0.0) tl[i] = true;
        else if (re > Lval && im == 0.0) tr[i] = true;
        else if (re < 0.0 && im > 0.0) bl[i] = true;
        else if (re < 0.0 && im == 0.0) brr[i] = true;
    }

    Eigen::VectorXcd zrect = Eigen::VectorXcd::Zero(n);
    for (int i = 0; i < 4; ++i) zrect(cnr[i] - 1) = rect[i];
    {
        const int len_ = c4 - c3 + 1;
        for (int j = 0; j < len_; ++j) zrect(c3 - 1 + j) = rect[2] + (rect[3] - rect[2]) * (static_cast<double>(j) / (len_ - 1));
    }
    {
        const int len_ = c2 - c1 + 1;
        for (int j = 0; j < len_; ++j) zrect(c1 - 1 + j) = rect[0] + (rect[1] - rect[0]) * (static_cast<double>(j) / (len_ - 1));
    }
    const double h = K / 20.0;
    {
        int cnt = 0;
        for (bool b : tl) if (b) ++cnt;
        int j = 1;
        for (int i = 0; i < n; ++i)
            if (tl[i]) {
                zrect(i) = std::complex<double>(0.0, Kp) - h * static_cast<double>(j);
                ++j;
            }
    }
    {
        int cnt = 0;
        for (bool b : tr) if (b) ++cnt;
        int j = 1;
        for (int i = 0; i < n; ++i)
            if (tr[i]) {
                zrect(i) = std::complex<double>(0.0, Kp) + h * static_cast<double>(cnt - j + 1);
                ++j;
            }
    }
    {
        int cnt = 0;
        for (bool b : bl) if (b) ++cnt;
        int j = 1;
        for (int i = 0; i < n; ++i)
            if (bl[i]) {
                zrect(i) = -h * static_cast<double>(cnt - j + 1);
                ++j;
            }
    }
    {
        int j = 1;
        for (int i = 0; i < n; ++i)
            if (brr[i]) {
                zrect(i) = h * static_cast<double>(j);
                ++j;
            }
    }

    Eigen::VectorXcd zn = zrect;
    std::vector<bool> done(n, false);
    for (int i = 0; i < 4; ++i) done[cnr[i] - 1] = true;
    const Eigen::VectorXcd cornerZ = Eigen::Map<const Eigen::VectorXcd>(rect.data(), 4);

    int iter = 0;
    const int maxiter = 50;
    while (iter < maxiter) {
        int remaining = 0;
        for (bool d : done)
            if (!d) ++remaining;
        if (remaining == 0) break;

        std::vector<int> active;
        for (int i = 0; i < n; ++i)
            if (!done[i]) active.push_back(i);
        const int mm = static_cast<int>(active.size());
        Eigen::VectorXcd znActive(mm);
        for (int j = 0; j < mm; ++j) znActive(j) = zn(active[j]);

        const R2Strip r2s = r2strip(znActive, cornerZ, Lval);
        Eigen::VectorXcd F(mm), step(mm);
        for (int j = 0; j < mm; ++j) {
            F(j) = zStrip(active[j]) - r2s.yp(j);
            step(j) = F(j) / r2s.yprime(j);
        }
        for (int j = 0; j < mm; ++j) {
            const int i = active[j];
            if (rSide[i] || lSide[i]) step(j) = std::complex<double>(0.0, step(j).imag());
            else step(j) = std::complex<double>(step(j).real(), 0.0);
        }
        Eigen::VectorXcd znew(mm);
        for (int j = 0; j < mm; ++j) znew(j) = znActive(j) + step(j);

        for (int j = 0; j < mm; ++j) {
            const int i = active[j];
            if (rSide[i] || lSide[i]) {
                const double x = std::min(std::max(znew(j).imag(), 0.0), Kp);
                znew(j) = std::complex<double>(znew(j).real(), x);
            }
        }
        for (int j = 0; j < mm; ++j) {
            const int i = active[j];
            if (tl[i] || bl[i]) {
                const double x = std::min(std::max(znew(j).real(), -K), -std::numeric_limits<double>::epsilon());
                znew(j) = std::complex<double>(x, znew(j).imag());
            }
        }
        for (int j = 0; j < mm; ++j) {
            const int i = active[j];
            if (tr[i] || brr[i]) {
                const double x = std::min(std::max(znew(j).real(), std::numeric_limits<double>::epsilon()), K);
                znew(j) = std::complex<double>(x, znew(j).imag());
            }
        }

        for (int j = 0; j < mm; ++j) zn(active[j]) = znew(j);
        for (int j = 0; j < mm; ++j)
            if (std::abs(F(j)) < tol) done[active[j]] = true;
        ++iter;
    }

    Eigen::VectorXcd zFinal(n);
    for (int i = 0; i < n; ++i) zFinal(renum[i]) = zn(i);

    Eigen::MatrixXd qdatOrig = qdat;
    for (int i = 0; i < n; ++i) {
        qdatOrig.col(renum[i]) = qdat.col(i);
        qdatOrig.col(n + 1 + renum[i]) = qdat.col(n + 1 + i);
    }

    return RParamResult{zFinal, c, Lval, qdatOrig};
}

}  // namespace sctoolbox
