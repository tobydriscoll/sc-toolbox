#include "sctoolbox/stparam.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "sctoolbox/nesolve.hpp"
#include "sctoolbox/scqdata.hpp"
#include "sctoolbox/stquad.hpp"
#include "sctoolbox/stquadh.hpp"

namespace sctoolbox {

namespace {

// Converts y (length n-1) into the n finite prevertices, matching
// stparam.m's/stpfun.m's shared formula:
// z(2:nb) = cumsum(exp(y(1:nb-1))); z(nb+1:n) = i+cumsum([y(nb);-exp(y(nb+1:n-1))]).
Eigen::VectorXcd yToZSt(const Eigen::VectorXd& y, int n, int nb) {
    Eigen::VectorXcd z = Eigen::VectorXcd::Zero(n);
    double cum = 0.0;
    for (int j = 0; j < nb - 1; ++j) {
        cum += std::exp(y(j));
        z(1 + j) = cum;
    }
    double cs = y(nb - 1);
    z(nb) = std::complex<double>(0.0, 1.0) + cs;
    for (int t = 1; t < n - nb; ++t) {
        cs -= std::exp(y(nb - 1 + t));
        z(nb + t) = std::complex<double>(0.0, 1.0) + cs;
    }
    return z;
}

// Port of stparam.m's nested stpfun residual function. left/right are
// 1-indexed finite-z positions (length n-1), already filtered/forced as
// in the outer stparam function. betaFull has length n+2 (the renumbered,
// un-reduced beta -- includes the two strip-end vertices' real turning
// angles, paired with the augmented zs array for stquad/stquadh).
Eigen::VectorXd stpfun(const Eigen::VectorXd& y, int n, int nb, const Eigen::VectorXd& betaFull,
                       const Eigen::VectorXcd& nmlen, const std::vector<int>& left, const std::vector<int>& right,
                       const std::vector<bool>& cmplx, const Eigen::MatrixXd& qdat) {
    const Eigen::VectorXcd z = yToZSt(y, n, nb);
    const int m = static_cast<int>(left.size());

    Eigen::VectorXcd zleft(m), zright(m), mid(m);
    for (int i = 0; i < m; ++i) {
        zleft(i) = z(left[i] - 1);
        zright(i) = z(right[i] - 1);
        mid(i) = (zleft(i) + zright(i)) / 2.0;
    }

    std::vector<bool> c2a = cmplx;
    c2a[1] = false;  // c2(2) = 0
    for (int i = 0; i < m; ++i) {
        if (c2a[i]) {
            const double diff = static_cast<double>(left[i] - nb);
            const double s = (diff > 0.0) ? 1.0 : (diff < 0.0 ? -1.0 : 0.0);
            mid(i) -= std::complex<double>(0.0, s / 2.0);
        }
    }

    const double inf = std::numeric_limits<double>::infinity();
    Eigen::VectorXcd zs(n + 2);
    zs(0) = std::complex<double>(-inf, 0.0);
    for (int j = 0; j < nb; ++j) zs(1 + j) = z(j);
    zs(nb + 1) = std::complex<double>(inf, 0.0);
    for (int j = 0; j < n - nb; ++j) zs(nb + 2 + j) = z(nb + j);

    std::vector<int> leftZs(m), rightZs(m);
    for (int i = 0; i < m; ++i) {
        leftZs[i] = left[i] + 1 + (left[i] > nb ? 1 : 0);
        rightZs[i] = right[i] + 1 + (right[i] > nb ? 1 : 0);
    }

    std::vector<bool> c2b = c2a;
    c2b[0] = true;  // c2(1) = 1

    std::vector<int> idsOnSide, idsAcross;
    for (int i = 0; i < m; ++i) (c2b[i] ? idsAcross : idsOnSide).push_back(i);

    Eigen::VectorXcd ints = Eigen::VectorXcd::Zero(m);

    if (!idsOnSide.empty()) {
        const int mm = static_cast<int>(idsOnSide.size());
        Eigen::VectorXcd zlv(mm), midv(mm), zrv(mm);
        std::vector<int> singL(mm), singR(mm);
        for (int j = 0; j < mm; ++j) {
            const int i = idsOnSide[j];
            zlv(j) = zleft(i);
            zrv(j) = zright(i);
            midv(j) = mid(i);
            singL[j] = leftZs[i];
            singR[j] = rightZs[i];
        }
        const Eigen::VectorXcd I1 = stquadh(zlv, midv, singL, zs, betaFull, qdat);
        const Eigen::VectorXcd I2 = stquadh(zrv, midv, singR, zs, betaFull, qdat);
        for (int j = 0; j < mm; ++j) ints(idsOnSide[j]) = I1(j) - I2(j);
    }

    if (!idsAcross.empty()) {
        const int mm = static_cast<int>(idsAcross.size());
        Eigen::VectorXcd zlv(mm), zrv(mm), z1v(mm), z2v(mm);
        std::vector<int> singL(mm), singR(mm), zeroSing(mm, 0);
        for (int j = 0; j < mm; ++j) {
            const int i = idsAcross[j];
            zlv(j) = zleft(i);
            zrv(j) = zright(i);
            z1v(j) = std::complex<double>(zleft(i).real(), 0.5);
            z2v(j) = std::complex<double>(zright(i).real(), 0.5);
            singL[j] = leftZs[i];
            singR[j] = rightZs[i];
        }
        const Eigen::VectorXcd I1 = stquad(zlv, z1v, singL, zs, betaFull, qdat);
        const Eigen::VectorXcd I2 = stquadh(z1v, z2v, zeroSing, zs, betaFull, qdat);
        const Eigen::VectorXcd I3 = stquad(zrv, z2v, singR, zs, betaFull, qdat);
        for (int j = 0; j < mm; ++j) ints(idsAcross[j]) = I1(j) + I2(j) - I3(j);
    }

    std::vector<int> idxF1, idxF2;
    for (int i = 0; i < m; ++i) (cmplx[i] ? idxF2 : idxF1).push_back(i);

    const double absval0 = std::abs(ints(idxF1[0]));
    const int n1 = static_cast<int>(idxF1.size());
    const int n2 = static_cast<int>(idxF2.size());

    Eigen::VectorXd rat1(n1 - 1);
    Eigen::VectorXcd rat2(n2);
    if (absval0 == 0.0) {
        rat1.setZero();
        rat2.setZero();
    } else {
        for (int j = 1; j < n1; ++j) rat1(j - 1) = std::abs(ints(idxF1[j])) / absval0;
        for (int j = 0; j < n2; ++j) rat2(j) = ints(idxF2[j]) / ints(idxF1[0]);
    }

    // cmplx2 = cmplx(2:end) (1-indexed) -- the per-pair complex flag for
    // entries 2..m (1-indexed), used to split nmlen (which already has its
    // own first entry dropped) into the F1/F2 groups by the SAME pattern.
    std::vector<bool> cmplx2(cmplx.begin() + 1, cmplx.end());
    std::vector<int> nmIdxF1, nmIdxF2;
    for (int i = 0; i < static_cast<int>(cmplx2.size()); ++i) (cmplx2[i] ? nmIdxF2 : nmIdxF1).push_back(i);

    Eigen::VectorXd F1(rat1.size());
    for (int j = 0; j < rat1.size(); ++j) F1(j) = std::log(rat1(j) / nmlen(nmIdxF1[j]).real());
    Eigen::VectorXcd F2(rat2.size());
    for (int j = 0; j < rat2.size(); ++j) F2(j) = std::log(rat2(j) / nmlen(nmIdxF2[j]));

    Eigen::VectorXd F(F1.size() + 2 * F2.size());
    int pos = 0;
    for (int j = 0; j < F1.size(); ++j) F(pos++) = F1(j);
    for (int j = 0; j < F2.size(); ++j) F(pos++) = F2(j).real();
    for (int j = 0; j < F2.size(); ++j) F(pos++) = F2(j).imag();
    return F;
}

}  // namespace

StParamResult stparam(const Eigen::VectorXcd& wIn, const Eigen::VectorXd& betaIn, std::array<int, 2> ends,
                      double tol, int method) {
    const int N = static_cast<int>(wIn.size());

    std::vector<int> renum(N);  // 0-indexed vertex ids, in renumbered order
    const int ends0 = ends[0] - 1;
    for (int i = 0; i < N; ++i) renum[i] = (ends0 + i) % N;

    Eigen::VectorXcd w(N);
    Eigen::VectorXd beta(N);
    for (int i = 0; i < N; ++i) {
        w(i) = wIn(renum[i]);
        beta(i) = betaIn(renum[i]);
    }

    int k0 = -1;
    const int ends2_0 = ends[1] - 1;
    for (int i = 0; i < N; ++i)
        if (renum[i] == ends2_0) {
            k0 = i;
            break;
        }
    const int k = k0 + 1;  // 1-indexed, matching MATLAB's k

    const int n = N - 2;
    const int nb = k - 2;

    const int nqpts = std::max(static_cast<int>(std::ceil(-std::log10(tol))), 4);
    const Eigen::MatrixXd qdatFull = scqdata(beta, nqpts);  // N+1 column-pairs

    std::vector<bool> atinf(N);
    for (int i = 0; i < N; ++i) atinf[i] = beta(i) <= -1.0;

    // w([1,k]) = []; atinf([1,k]) = []; (1-indexed positions 1 and k)
    Eigen::VectorXcd wRed(n);
    std::vector<bool> atinfRed(n);
    {
        int p = 0;
        for (int i = 0; i < N; ++i) {
            if (i == 0 || i == k - 1) continue;
            wRed(p) = w(i);
            atinfRed[p] = atinf[i];
            ++p;
        }
    }

    // Initial guess z0 (length n), matching the z0=[] branch of stparam.m.
    Eigen::VectorXcd z0(n);
    bool anyAtinf = false;
    for (bool b : atinfRed)
        if (b) anyAtinf = true;

    if (anyAtinf) {
        const double scale = (std::abs(wRed(nb - 1) - wRed(0)) + std::abs(wRed(n - 1) - wRed(nb))) / 2.0;
        for (int j = 0; j < nb; ++j) z0(j) = scale * j / (nb - 1);
        for (int j = 0; j < n - nb; ++j) z0(nb + j) = std::complex<double>(0.0, 1.0) + scale * (n - nb - 1 - j) / (n - nb - 1);
    } else {
        const double scale1 = (std::abs(wRed(n - 1) - wRed(0)) + std::abs(wRed(nb - 1) - wRed(nb))) / 2.0;
        Eigen::VectorXd z0r(n);
        z0r(0) = 0.0;
        for (int j = 1; j < nb; ++j) z0r(j) = z0r(j - 1) + std::abs(wRed(j) - wRed(j - 1)) / scale1;
        if (nb + 1 == n) {
            z0r(n - 1) = (z0r(0) + z0r(nb - 1)) / 2.0;
        } else {
            z0r(n - 1) = 0.0;
            for (int j = n - 2; j >= nb; --j) z0r(j) = z0r(j + 1) + std::abs(wRed(j) - wRed(j + 1)) / scale1;
        }
        const double scale2 = std::sqrt(z0r(nb - 1) / z0r(nb));
        for (int j = 0; j < nb; ++j) z0(j) = z0r(j) / scale2;
        for (int j = nb; j < n; ++j) z0(j) = std::complex<double>(0.0, 1.0) + z0r(j) * scale2;
    }

    Eigen::VectorXd y0(n - 1);
    for (int j = 0; j < nb - 1; ++j) y0(j) = std::log((z0(j + 1) - z0(j)).real());
    y0(nb - 1) = z0(nb).real();
    for (int j = 0; j < n - nb - 1; ++j) y0(nb + j) = std::log(-(z0(nb + 1 + j) - z0(nb + j)).real());

    // left/right (1-indexed finite-z positions), then deletions.
    std::vector<int> leftArr, rightArr;
    leftArr.push_back(1);
    for (int i = 2; i <= n; ++i) leftArr.push_back(i - 1);
    rightArr.push_back(n);
    for (int i = 2; i <= n; ++i) rightArr.push_back(i);

    std::vector<int> delLeftPos, delRightPos;
    for (int p = 1; p <= n; ++p)
        if (atinfRed[p - 1]) delLeftPos.push_back(p + 1);
    delLeftPos.push_back(nb + 1);
    for (int p = 1; p <= n; ++p)
        if (atinfRed[p - 1]) delRightPos.push_back(p);
    delRightPos.push_back(nb + 1);

    auto eraseAt = [](std::vector<int>& arr, std::vector<int> positions) {
        std::sort(positions.rbegin(), positions.rend());
        for (int pos : positions) arr.erase(arr.begin() + (pos - 1));
    };
    eraseAt(leftArr, delLeftPos);
    eraseAt(rightArr, delRightPos);

    const int m = static_cast<int>(leftArr.size());
    std::vector<bool> cmplx(m);
    for (int i = 0; i < m; ++i) cmplx[i] = (rightArr[i] - leftArr[i] == 2);
    cmplx[0] = false;
    cmplx[1] = true;

    Eigen::VectorXcd nmlenC(m);
    for (int i = 0; i < m; ++i) {
        nmlenC(i) = (wRed(rightArr[i] - 1) - wRed(leftArr[i] - 1)) / (wRed(n - 1) - wRed(0));
        if (!cmplx[i]) nmlenC(i) = std::abs(nmlenC(i));
    }
    // nmlen stays complex (the cmplx-true entries keep their genuine
    // imaginary part, needed by F2's complex division/log in stpfun); only
    // the cmplx-false entries are real-valued (abs'd above).
    const Eigen::VectorXcd nmlen = nmlenC.tail(m - 1);

    const Fvec fvec = [&](const Eigen::VectorXd& y) {
        return stpfun(y, n, nb, beta, nmlen, leftArr, rightArr, cmplx, qdatFull);
    };

    Eigen::VectorXd details = Eigen::VectorXd::Zero(16);
    details(0) = 0.0;
    details(1) = static_cast<double>(method);
    details(5) = 100.0 * (n - 1);
    details(7) = tol;
    details(8) = tol / 10.0;
    details(11) = nqpts;

    const NesolveResult r = nesolve(fvec, y0, details);
    const Eigen::VectorXcd zFinite = yToZSt(r.xf, n, nb);

    const double inf = std::numeric_limits<double>::infinity();
    Eigen::VectorXcd zAug(N);
    zAug(0) = std::complex<double>(-inf, 0.0);
    for (int j = 0; j < nb; ++j) zAug(1 + j) = zFinite(j);
    zAug(nb + 1) = std::complex<double>(inf, 0.0);
    for (int j = 0; j < n - nb; ++j) zAug(nb + 2 + j) = zFinite(nb + j);

    const std::complex<double> mid = (zAug(1) + zAug(2)) / 2.0;
    const std::vector<int> singAt3{3}, singAt2{2};
    const Eigen::VectorXcd g1 = stquad(zAug.segment(2, 1), Eigen::VectorXcd::Constant(1, mid), singAt3, zAug, beta, qdatFull);
    const Eigen::VectorXcd g2 = stquad(zAug.segment(1, 1), Eigen::VectorXcd::Constant(1, mid), singAt2, zAug, beta, qdatFull);
    const std::complex<double> g = g1(0) - g2(0);
    const std::complex<double> c = (wRed(0) - wRed(1)) / g;

    // Undo renumbering.
    Eigen::VectorXcd zOrig(N);
    for (int i = 0; i < N; ++i) zOrig(renum[i]) = zAug(i);

    Eigen::MatrixXd qdatOrig = qdatFull;
    for (int i = 0; i < N; ++i) {
        qdatOrig.col(renum[i]) = qdatFull.col(i);
        qdatOrig.col(N + 1 + renum[i]) = qdatFull.col(N + 1 + i);
    }

    return StParamResult{zOrig, c, qdatOrig};
}

}  // namespace sctoolbox
