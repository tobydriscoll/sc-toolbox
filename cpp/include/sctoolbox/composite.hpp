#pragma once
#include <Eigen/Dense>
#include <functional>
#include <string>
#include <utility>
#include <vector>

#include "sctoolbox/moebius.hpp"

namespace sctoolbox {

// Port of @composite: the map obtained by applying a list of maps in turn,
// f = fN(...f2(f1(z))...).
//
// MATLAB's composite.m relies on duck typing -- it accepts any of moebius,
// the six SC map classes, scmapinv, and inline functions, and dispatches
// through feval. C++ needs explicit type erasure instead, so each member is
// stored as a forward callable plus (when the member is invertible) its
// inverse. That pairing is what makes Composite::inverse() work: MATLAB's
// inv() reverses the list and calls inv() on each member, which for an SC
// map yields an scmapinv whose eval is the map's evalinv -- exactly the
// forward/inverse swap performed here.
//
// Build members with the factories below rather than by hand:
//
//   Composite f{Composite::member(dmap), Composite::member(mob)};
//   Eigen::VectorXcd w = f.eval(z);
//   Composite fi = f.inverse();
class Composite {
public:
    using Fn = std::function<Eigen::VectorXcd(const Eigen::VectorXcd&)>;

    struct Member {
        Fn forward;
        Fn inverse;  // empty for members that cannot be inverted
        std::string name;
    };

    Composite() = default;
    Composite(std::initializer_list<Member> members) : maps_(members) {}
    explicit Composite(std::vector<Member> members) : maps_(std::move(members)) {}

    // MATLAB's constructor flattens a nested composite into its members.
    Composite& append(Member m);
    Composite& append(const Composite& other);

    // Port of @composite/eval.m.
    Eigen::VectorXcd eval(const Eigen::VectorXcd& z) const;

    // Port of @composite/inv.m. Throws std::runtime_error if any member has
    // no inverse, mirroring MATLAB's "Can't invert INLINE maps." error.
    Composite inverse() const;

    // Port of @composite/members.m.
    const std::vector<Member>& members() const { return maps_; }
    int length() const { return static_cast<int>(maps_.size()); }

    // Any SC map class with eval/evalinv (DiskMap, HplMap, ExterMap,
    // StripMap, RectMap, CrDiskMap). The map is copied into the member, so
    // the composite stays valid independently of the caller's object.
    template <class Map>
    static Member member(const Map& m, std::string name = "scmap") {
        Member out;
        out.forward = [m](const Eigen::VectorXcd& z) { return m.eval(z); };
        out.inverse = [m](const Eigen::VectorXcd& w) { return m.evalinv(w).zp; };
        out.name = std::move(name);
        return out;
    }

    // The scmapinv of an SC map: forward and inverse exchanged.
    template <class Map>
    static Member inverseMember(const Map& m, std::string name = "scmapinv") {
        Member out = member(m, std::move(name));
        std::swap(out.forward, out.inverse);
        return out;
    }

    static Member member(const Moebius& m, std::string name = "moebius");

    // An arbitrary function, standing in for MATLAB's inline members: usable
    // in a composite but not invertible.
    static Member function(Fn f, std::string name = "function");

private:
    std::vector<Member> maps_;
};

}  // namespace sctoolbox
