#include "sctoolbox/composite.hpp"

#include <stdexcept>

namespace sctoolbox {

Composite& Composite::append(Member m) {
    maps_.push_back(std::move(m));
    return *this;
}

Composite& Composite::append(const Composite& other) {
    maps_.insert(maps_.end(), other.maps_.begin(), other.maps_.end());
    return *this;
}

Eigen::VectorXcd Composite::eval(const Eigen::VectorXcd& z) const {
    Eigen::VectorXcd w = z;
    for (const auto& m : maps_) w = m.forward(w);
    return w;
}

Composite Composite::inverse() const {
    std::vector<Member> list;
    list.reserve(maps_.size());
    for (auto it = maps_.rbegin(); it != maps_.rend(); ++it) {
        if (!it->inverse) {
            throw std::runtime_error("Composite::inverse: member '" + it->name +
                                     "' has no inverse.");
        }
        Member m = *it;
        std::swap(m.forward, m.inverse);
        list.push_back(std::move(m));
    }
    return Composite(std::move(list));
}

Composite::Member Composite::member(const Moebius& m, std::string name) {
    Member out;
    out.forward = [m](const Eigen::VectorXcd& z) { return m.eval(z); };
    const Moebius mi = m.inverse();
    out.inverse = [mi](const Eigen::VectorXcd& w) { return mi.eval(w); };
    out.name = std::move(name);
    return out;
}

Composite::Member Composite::function(Fn f, std::string name) {
    Member out;
    out.forward = std::move(f);
    out.name = std::move(name);
    return out;
}

}  // namespace sctoolbox
