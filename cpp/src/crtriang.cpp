#include "sctoolbox/crtriang.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include "sctoolbox/crpsdist.hpp"
#include "sctoolbox/isinpoly.hpp"
#include "sctoolbox/scangle.hpp"

namespace sctoolbox {

namespace {
// Finds the 1-indexed column of the first column in M(:, 1:lastCol) that is
// entirely zero, matching MATLAB's `find(any(M==0), 1)`. Returns -1 if none.
int firstZeroColumn(const Eigen::MatrixXi& M, int lastCol) {
    for (int j = 0; j < lastCol; ++j)
        if ((M.col(j).array() == 0).any()) return j + 1;
    return -1;
}
}  // namespace

CrTriangulation crtriang(const Eigen::VectorXcd& Win) {
    const int N = static_cast<int>(Win.size());
    Eigen::MatrixXi edge = Eigen::MatrixXi::Zero(2, 2 * N - 3);
    Eigen::MatrixXi triedge = Eigen::MatrixXi::Zero(3, N - 2);
    Eigen::MatrixXi edgetri = Eigen::MatrixXi::Zero(2, 2 * N - 3);

    for (int i = 1; i <= N; ++i) {
        edge(0, N - 3 + i - 1) = i;
        edge(1, N - 3 + i - 1) = (i % N) + 1;
    }

    std::vector<std::vector<int>> stack;
    {
        std::vector<int> all(N);
        for (int i = 0; i < N; ++i) all[i] = i + 1;
        stack.push_back(all);
    }

    const double epsD = std::numeric_limits<double>::epsilon();

    while (!stack.empty()) {
        const std::vector<int> idx = stack.back();
        stack.pop_back();
        const int n = static_cast<int>(idx.size());

        Eigen::VectorXcd w(n);
        for (int i = 0; i < n; ++i) w(i) = Win(idx[i] - 1);

        if (n == 3) {
            const int tnum = firstZeroColumn(triedge, N - 2);
            int e[3][2];
            const int ev[6] = {idx[0], idx[1], idx[1], idx[2], idx[2], idx[0]};
            for (int j = 0; j < 3; ++j) {
                int a = ev[2 * j], b = ev[2 * j + 1];
                if (a > b) std::swap(a, b);
                e[j][0] = a;
                e[j][1] = b;
            }
            for (int j = 0; j < 3; ++j) {
                const int de = e[j][1] - e[j][0];
                int enumv;
                if (de == 1) {
                    enumv = N - 3 + e[j][0];
                } else if (de == N - 1) {
                    enumv = 2 * N - 3;
                } else {
                    enumv = -1;
                    for (int c = 0; c < N - 3; ++c) {
                        if (edge(0, c) == e[j][0] && edge(1, c) == e[j][1]) {
                            enumv = c + 1;
                            break;
                        }
                    }
                }
                if (enumv < 1) {
                    // Should be unreachable now that edges are stored sorted
                    // (see the diagonal-storage note below). Guard defensively
                    // so a lookup miss surfaces as a catchable error rather than
                    // an out-of-bounds write that corrupts the heap.
                    throw std::runtime_error(
                        "crtriang: could not resolve a triangle edge during "
                        "polygon triangulation");
                }
                triedge(j, tnum - 1) = enumv;
                if (edgetri(0, enumv - 1) == 0)
                    edgetri(0, enumv - 1) = tnum;
                else
                    edgetri(1, enumv - 1) = tnum;
            }
        } else {
            const Eigen::VectorXd beta = scangle(w);
            int j = 0;
            {
                double best = beta(0);
                for (int i = 1; i < n; ++i)
                    if (beta(i) < best) {
                        best = beta(i);
                        j = i;
                    }
            }
            const int jm1 = (j + n - 1) % n;
            const int jp1 = (j + 1) % n;

            std::vector<bool> t(n, false);
            t[jm1] = t[j] = t[jp1] = true;

            Eigen::VectorXcd triangle(3);
            triangle(0) = w(jm1);
            triangle(1) = w(j);
            triangle(2) = w(jp1);

            std::vector<int> others;
            for (int i = 0; i < n; ++i)
                if (!t[i]) others.push_back(i);
            const int no = static_cast<int>(others.size());

            std::vector<bool> inside(n, false);
            if (no > 0) {
                Eigen::VectorXcd wOthers(no);
                for (int i = 0; i < no; ++i) wOthers(i) = w(others[i]);
                const Eigen::VectorXd ip = isinpoly(wOthers, triangle);
                for (int i = 0; i < no; ++i) inside[others[i]] = (ip(i) != 0.0);
            }

            for (int k = 0; k < 3; ++k) {
                const std::complex<double> segA = triangle((k) % 3);
                const std::complex<double> segB = triangle((k + 1) % 3);
                if (no > 0) {
                    Eigen::VectorXcd wOthers(no);
                    for (int i = 0; i < no; ++i) wOthers(i) = w(others[i]);
                    const Eigen::VectorXd d = crpsdist({segA, segB}, wOthers);
                    for (int i = 0; i < no; ++i) {
                        if (d(i) < epsD) {
                            const int p = others[i];
                            inside[p] = false;
                            double minDistToTri = std::numeric_limits<double>::infinity();
                            for (int q = 0; q < 3; ++q) minDistToTri = std::min(minDistToTri, std::abs(w(p) - triangle(q)));
                            if (minDistToTri > epsD) {
                                const int pPrev = (p - 1 + n) % n;
                                const std::complex<double> s =
                                    std::exp(std::complex<double>(0.0, std::arg(w(pPrev) - w(p)) - M_PI * (beta(p) + 1.0) / 2.0));
                                Eigen::VectorXcd testPt(1);
                                testPt(0) = w(p) + 1e-13 * std::abs(w(p)) * s;
                                const Eigen::VectorXd r = isinpoly(testPt, triangle);
                                if (r(0) != 0.0) inside[p] = true;
                            }
                        }
                    }
                }
            }

            std::vector<int> insideIdx;
            for (int i = 0; i < n; ++i)
                if (inside[i]) insideIdx.push_back(i);

            int e1, e2;  // 0-indexed positions within the sub-polygon
            if (insideIdx.empty()) {
                e1 = jm1;
                e2 = jp1;
                if (e1 > e2) std::swap(e1, e2);
            } else {
                int chosen;
                if (insideIdx.size() > 1) {
                    const int m = static_cast<int>(insideIdx.size());
                    std::vector<int> vis;
                    for (int ii = 0; ii < m; ++ii) {
                        const int p = insideIdx[ii];
                        const int fwd = (p + 1) % n;
                        double ang1 = std::arg((w(j) - w(p)) / (w(fwd) - w(p)));
                        ang1 = std::fmod(ang1 / M_PI + 2.0, 2.0);
                        double ang2 = std::arg((w(p) - w(j)) / (w(jp1) - w(j)));
                        ang2 = std::fmod(ang2 / M_PI + 2.0, 2.0);
                        if (ang1 < beta(p) + 1.0 && ang2 < beta(j) + 1.0) vis.push_back(ii);
                    }
                    const std::complex<double> dw1 = w(jp1) - w(j);
                    const std::complex<double> dw2 = w(j) - w(jm1);
                    const double theta = std::arg(dw1) + std::arg(dw2 / dw1) / 2.0;
                    const std::complex<double> rot = std::exp(std::complex<double>(0.0, theta));
                    double bestD = std::numeric_limits<double>::infinity();
                    int bestVis = vis.empty() ? 0 : vis[0];
                    for (int v : vis) {
                        const std::complex<double> wvis = w(insideIdx[v]) - w(j);
                        const double D = std::abs(wvis - (wvis * std::conj(rot)).real() * rot);
                        if (D < bestD) {
                            bestD = D;
                            bestVis = v;
                        }
                    }
                    chosen = insideIdx[bestVis];
                } else {
                    chosen = insideIdx[0];
                }
                e1 = chosen;
                e2 = j;
                if (e1 > e2) std::swap(e1, e2);
            }

            const int enumv = firstZeroColumn(edge, 2 * N - 3);
            // Store the diagonal's endpoints sorted by *original* vertex index.
            // MATLAB represents each sub-polygon as a boolean mask and re-reads
            // it with find(), so its vertex list is always in increasing global
            // order and `edge(:,enum)=idx(e)` is implicitly sorted. This port
            // keeps a rotated vertex list (i2 wraps around, see below), so e1<e2
            // by *position* does not imply idx[e1]<idx[e2]; sorting here restores
            // the invariant the base-case edge lookup (line ~72) relies on.
            // Without this, polygons that subdivide past a single split (n>=5)
            // record an unsorted edge, the lookup misses it, and enumv stays -1.
            int va = idx[e1], vb = idx[e2];
            if (va > vb) std::swap(va, vb);
            edge(0, enumv - 1) = va;
            edge(1, enumv - 1) = vb;

            std::vector<int> i1, i2;
            for (int p = e1; p <= e2; ++p) i1.push_back(idx[p]);
            for (int p = e2; p < n; ++p) i2.push_back(idx[p]);
            for (int p = 0; p <= e1; ++p) i2.push_back(idx[p]);

            stack.push_back(i1);
            stack.push_back(i2);
        }
    }

    return CrTriangulation{edge, triedge, edgetri};
}

}  // namespace sctoolbox
