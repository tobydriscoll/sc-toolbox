#include "sctoolbox/scfix.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>

#include "sctoolbox/scaddvtx.hpp"
#include "sctoolbox/scangle.hpp"

namespace sctoolbox {

namespace {

const double kEps = std::numeric_limits<double>::epsilon();

bool isInf(const std::complex<double>& z) { return std::isinf(z.real()) || std::isinf(z.imag()); }

// All permutation/index bookkeeping below mirrors scfix.m's 1-indexed
// arithmetic directly (P[i-1] holds the 1-indexed source position for new
// position i) so that the control flow can be checked line-by-line against
// the MATLAB source; only the final array reads/writes drop to 0-indexed.

std::vector<int> rotateFrom(int p, int n) {
    std::vector<int> r(n);
    int idx = 0;
    for (int i = p; i <= n; ++i) r[idx++] = i;
    for (int i = 1; i < p; ++i) r[idx++] = i;
    return r;
}

Eigen::VectorXcd permute1(const Eigen::VectorXcd& v, const std::vector<int>& p) {
    const int n = static_cast<int>(p.size());
    Eigen::VectorXcd out(n);
    for (int i = 0; i < n; ++i) out(i) = v(p[i] - 1);
    return out;
}

Eigen::VectorXd permute1d(const Eigen::VectorXd& v, const std::vector<int>& p) {
    const int n = static_cast<int>(p.size());
    Eigen::VectorXd out(n);
    for (int i = 0; i < n; ++i) out(i) = v(p[i] - 1);
    return out;
}

std::vector<int> permute1i(const std::vector<int>& v, const std::vector<int>& p) {
    const int n = static_cast<int>(p.size());
    std::vector<int> out(n);
    for (int i = 0; i < n; ++i) out[i] = v[p[i] - 1];
    return out;
}

int find1(const std::vector<int>& v, int val) {
    for (size_t i = 0; i < v.size(); ++i)
        if (v[i] == val) return static_cast<int>(i) + 1;
    throw std::runtime_error("scfix: value not found in renum array");
}

}  // namespace

ScfixResult scfix(const std::string& type, Eigen::VectorXcd w, Eigen::VectorXd beta, std::vector<int> aux) {
    int n = static_cast<int>(w.size());
    std::vector<int> renum(n);
    for (int i = 0; i < n; ++i) renum[i] = i + 1;

    // Orientation convention: reverse traversal order (fixing vertex 1) if
    // the supplied beta sum has the wrong sign.
    const double sumb = -2.0 + (type == "de" ? 4.0 : 0.0);
    if (std::abs(beta.sum() + sumb) < 1e-9) {
        std::vector<int> idxRev(n);
        idxRev[0] = 1;
        for (int k = 2; k <= n; ++k) idxRev[k - 1] = n - k + 2;
        w = permute1(w, idxRev);
        beta = scangle(w);
        renum = permute1i(renum, idxRev);
        if (!aux.empty()) {
            std::vector<int> newAux(aux.size());
            for (size_t i = 0; i < aux.size(); ++i) newAux[i] = idxRev[aux[i] - 1];
            aux = newAux;
        }
    }

    if (type == "hp" || type == "d") {
        std::vector<int> shift(n);
        for (int i = 1; i < n; ++i) shift[i - 1] = i + 1;
        shift[n - 1] = 1;

        auto needsFix = [&]() {
            const bool infBad = isInf(w(0)) || isInf(w(1)) || isInf(w(n - 2));
            const bool betaBad = std::abs(beta(n - 1)) < kEps || std::abs(beta(n - 1) - 1.0) < kEps;
            return infBad || betaBad;
        };
        while (needsFix()) {
            renum = permute1i(renum, shift);
            w = permute1(w, shift);
            beta = permute1d(beta, shift);
            if (renum[0] == 1) {  // tried all orderings
                bool allDegenerate = true;
                for (int i = 0; i < n; ++i) {
                    if (!(std::abs(beta(i) - 1.0) < kEps || std::abs(beta(i)) < kEps)) {
                        allDegenerate = false;
                        break;
                    }
                }
                if (allDegenerate) throw std::runtime_error("Polygon has empty interior!");
                while (std::abs(beta(n - 1)) < kEps || std::abs(beta(n - 1) - 1.0) < kEps) {
                    w = permute1(w, shift);
                    beta = permute1d(beta, shift);
                }
                if (isInf(w(0)) || isInf(w(1))) {
                    auto r = scaddvtx(w, beta, 0);
                    w = r.w;
                    beta = r.beta;
                    n += 1;
                }
                if (isInf(w(n - 2))) {
                    auto r = scaddvtx(w, beta, n - 2);
                    w = r.w;
                    beta = r.beta;
                    n += 1;
                }
                renum.resize(n);
                for (int i = 0; i < n; ++i) renum[i] = i + 1;
                break;
            }
        }
    } else if (type == "de") {
        std::vector<int> shift(n);
        for (int i = 1; i < n; ++i) shift[i - 1] = i + 1;
        shift[n - 1] = 1;

        while ((std::abs(beta(n - 1)) < kEps || std::abs(beta(n - 1) - 1.0) < kEps) && n > 2) {
            renum = permute1i(renum, shift);
            w = permute1(w, shift);
            beta = permute1d(beta, shift);
            if (renum[0] == 1) {
                std::vector<std::complex<double>> wkeep;
                std::vector<double> bkeep;
                for (int i = 0; i < n; ++i) {
                    if (std::abs(beta(i)) >= kEps) {
                        wkeep.push_back(w(i));
                        bkeep.push_back(beta(i));
                    }
                }
                const int nkeep = static_cast<int>(wkeep.size());
                w.resize(nkeep);
                beta.resize(nkeep);
                for (int i = 0; i < nkeep; ++i) {
                    w(i) = wkeep[i];
                    beta(i) = bkeep[i];
                }
                renum = {1, 2};
                n = 2;
                break;
            }
        }
    } else if (type == "st") {
        if (aux.empty())
            throw std::runtime_error(
                "scfix: interactive strip-end selection not supported; aux (ends) must be provided");
        std::vector<int> ends = aux;

        auto perm = rotateFrom(ends[0], n);
        renum = perm;
        w = permute1(w, perm);
        beta = permute1d(beta, perm);
        int k = find1(renum, ends[1]);

        if (k < 4) {
            if (k < n - 1) {
                auto perm2 = rotateFrom(k, n);
                renum = perm2;
                w = permute1(w, perm2);
                beta = permute1d(beta, perm2);
                k = find1(renum, 1);
            } else {
                for (int j = 1; j <= 4 - k; ++j) {
                    auto r = scaddvtx(w, beta, j - 1);
                    w = r.w;
                    beta = r.beta;
                    n += 1;
                    k += 1;
                }
            }
        }

        if (k == n) {
            auto r = scaddvtx(w, beta, n - 1);
            w = r.w;
            beta = r.beta;
            n += 1;
        }

        if (isInf(w(1))) {
            for (int j = 1; j <= 2; ++j) {
                auto r = scaddvtx(w, beta, j - 1);
                w = r.w;
                beta = r.beta;
                n += 1;
                k += 1;
            }
        } else if (isInf(w(2))) {
            auto r = scaddvtx(w, beta, 1);
            w = r.w;
            beta = r.beta;
            n += 1;
            k += 1;
        } else if (isInf(w(n - 1))) {
            auto r = scaddvtx(w, beta, n - 1);
            w = r.w;
            beta = r.beta;
            n += 1;
        }

        aux = {1, k};
    } else if (type == "r") {
        if (aux.empty()) throw std::runtime_error("scfix: corner positions (aux) must be provided for type 'r'");
        std::vector<int> corner = aux;

        auto perm = rotateFrom(corner[0], n);
        renum = perm;
        w = permute1(w, perm);
        beta = permute1d(beta, perm);

        auto newPos = [](int c, int c1, int nn) { return ((c - c1) % nn + nn) % nn + 1; };
        const int c1 = corner[0];
        for (auto& c : corner) c = newPos(c, c1, n);

        if (std::abs(beta(n - 1)) < kEps || std::abs(beta(n - 1) - 1.0) < kEps) {
            const int prevLabel = corner[2] - 1;  // beta(corner(3)-1), 1-indexed label
            const bool cond =
                !(std::abs(beta(prevLabel - 1)) < kEps || std::abs(beta(prevLabel - 1) - 1.0) < kEps) &&
                !isInf(w(corner[2] - 1));
            if (cond) {
                auto perm2 = rotateFrom(corner[2], n);
                renum = perm2;
                w = permute1(w, perm2);
                beta = permute1d(beta, perm2);
                const int c3 = corner[2];
                std::vector<int> newCorner(4);
                for (int i = 0; i < 4; ++i) newCorner[i] = newPos(corner[i], c3, n);
                std::sort(newCorner.begin(), newCorner.end());
                corner = newCorner;
            } else {
                throw std::runtime_error("Collinear sides make posing problem impossible");
            }
        }

        if (isInf(w(1))) {
            auto r = scaddvtx(w, beta, 0);
            w = r.w;
            beta = r.beta;
            n += 1;
            for (int i = 1; i < 4; ++i) corner[i] += 1;
        }

        aux = corner;
    }

    return {w, beta, aux};
}

}  // namespace sctoolbox
