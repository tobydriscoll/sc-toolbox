#include "sctoolbox/rectmap_internal.hpp"

#include <vector>

namespace sctoolbox {

std::pair<int, int> rectStripEnds(const Eigen::VectorXcd& z) {
    const int n = static_cast<int>(z.size());
    int e0 = -1, e1 = -1;
    for (int i = 0; i < n; ++i) {
        const double d = z((i + 1) % n).imag() - z(i).imag();
        if (d != 0.0) {
            if (e0 < 0) e0 = i;
            else e1 = i;
        }
    }
    return {e0, e1};
}

Eigen::MatrixXd rectAugQdat(const Eigen::MatrixXd& qdat, int n, int e0, int e1) {
    // Matches MATLAB's `idx = [1:ends(1) n+1 ends(1)+1:ends(2) n+1 ends(2)+1:n n+1]`
    // (length n+3, with a *third* "n+1" filler that looks redundant but is not:
    // stmap/stquad index qdat assuming (n_aug+1)-sized column blocks -- the same
    // "neutral column" convention scqdata uses -- where n_aug = n+2 is the size of
    // the augmented zs/ws/bs arrays. So the augmented qdat needs n_aug+1 = n+3
    // node columns (and n+3 weight columns), not n+2.
    const int nqpts = static_cast<int>(qdat.rows());
    std::vector<int> idx;
    idx.reserve(n + 3);
    for (int i = 0; i <= e0; ++i) idx.push_back(i);
    idx.push_back(n);
    for (int i = e0 + 1; i <= e1; ++i) idx.push_back(i);
    idx.push_back(n);
    for (int i = e1 + 1; i < n; ++i) idx.push_back(i);
    idx.push_back(n);

    const int m = static_cast<int>(idx.size());  // n + 3
    Eigen::MatrixXd aug(nqpts, 2 * m);
    for (int k = 0; k < m; ++k) {
        aug.col(k) = qdat.col(idx[k]);
        aug.col(m + k) = qdat.col(n + 1 + idx[k]);
    }
    return aug;
}

}  // namespace sctoolbox
