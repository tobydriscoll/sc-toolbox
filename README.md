sc-toolbox
==========

A **C++17 port** of the Schwarz–Christoffel Toolbox for conformal mapping — numerical routines for computing Schwarz–Christoffel conformal maps onto regions bounded by polygons in the complex plane.

This is a fork of [tobydriscoll/sc-toolbox](https://github.com/tobydriscoll/sc-toolbox), the original MATLAB toolbox by Toby Driscoll. The MATLAB sources are still here and unchanged in behavior: they are the reference implementation, and they generate the golden values the C++ tests are checked against.

For the mathematics, see *Schwarz–Christoffel Mapping* by Driscoll and Trefethen. For a user's guide to the original toolbox (concepts, options, and map types all carry over), visit <https://tobydriscoll.net/project/sc-toolbox/>. The MATLAB version is also on the [File Exchange](https://www.mathworks.com/matlabcentral/fileexchange/1316-schwarz-christoffel-toolbox), where you can try it online.

## What's ported

Every map type in the toolbox, function by function, verified against MATLAB-generated golden values:

| C++ class | Domain | Target |
|---|---|---|
| `DiskMap` | unit disk | polygon interior |
| `HplMap` | upper half-plane | polygon interior |
| `ExterMap` | disk exterior | polygon exterior |
| `StripMap` | infinite strip | polygon interior |
| `RectMap` | rectangle | generalized quadrilateral |
| `CrDiskMap` | disk (cross-ratio formulation) | polygon interior |
| `AnnulusMap` | annulus | doubly connected polygonal region |

Plus `Polygon`, `Moebius`, `Composite`, the Gauss–Jacobi quadrature machinery, the `nesolve` Newton solver, and the elliptic functions the rectangle map needs. The port matches MATLAB to within `5e-11` absolute on all golden cases.

## Build

Prerequisites — [Eigen 3](https://eigen.tuxfamily.org), CMake 3.16+, and (for the tests only) [Catch2 3](https://github.com/catchorg/Catch2). On macOS:

```sh
brew install cmake eigen catch2
```

```sh
cmake -S cpp -B cpp/build -DCMAKE_PREFIX_PATH=/opt/homebrew
cmake --build cpp/build
./cpp/build/sctoolbox_tests     # or: ctest --test-dir cpp/build
```

MATLAB is **not** needed to build or test the C++ library — the golden values are committed as text. If Eigen is not installed, the build fetches a pinned copy automatically.

### Using it in your project

```cmake
add_subdirectory(path/to/sc-toolbox/cpp sctoolbox)
target_link_libraries(your_target PRIVATE sctoolbox)
```

The test suite is skipped automatically when the project is consumed this way, so Catch2 is not required. Force it either way with `-DSCTOOLBOX_BUILD_TESTS=ON/OFF`.

## Basic usage

```cpp
#include <iostream>
#include "sctoolbox/diskmap.hpp"

int main() {
    // The unit square, counterclockwise. Interior angles are computed from
    // the geometry; pass them explicitly for unbounded polygons.
    Eigen::VectorXcd v(4);
    v << std::complex<double>(0, 0), std::complex<double>(1, 0),
         std::complex<double>(1, 1), std::complex<double>(0, 1);
    const sctoolbox::Polygon square(v);

    // The parameter problem is solved in the constructor.
    const sctoolbox::DiskMap map(square);
    std::cout << "accuracy: " << map.accuracy() << "\n";

    Eigen::VectorXcd zp(3);
    zp << std::complex<double>(0.0, 0.0), std::complex<double>(0.5, 0.25),
          std::complex<double>(-0.3, 0.6);

    const Eigen::VectorXcd wp   = map.eval(zp);         // disk -> polygon
    const Eigen::VectorXcd back = map.evalinv(wp).zp;   // polygon -> disk
    const Eigen::VectorXcd dp   = map.evaldiff(zp);     // derivative

    for (int i = 0; i < zp.size(); ++i)
        std::cout << zp(i) << " -> " << wp(i) << " -> " << back(i)
                  << "   f'=" << dp(i) << "\n";
}
```

```
accuracy: 1.57996e-12
(0,0) -> (0.5,0.5) -> (2.07122e-12,2.07139e-12)   f'=(-0.38138,0.38138)
(0.5,0.25) -> (0.213966,0.592437) -> (0.5,0.25)   f'=(-0.392182,0.358023)
(-0.3,0.6) -> (0.392702,0.157038) -> (-0.3,0.6)   f'=(-0.399835,0.332987)
```

Every map class exposes the same four methods as its MATLAB counterpart: `eval`, `evalinv`, `evaldiff`, and `accuracy` (a self-consistency estimate, memoized after the first call). Points outside the domain come back as NaN from `eval`, matching MATLAB.

### Things to know

- **Constructors do the work.** Building a map solves the nonlinear parameter problem, which is the expensive step. Construct once and reuse; pass a tolerance as the last constructor argument (default `1e-8`).
- **`evalinv` returns a struct**, not a bare vector: `.zp` holds the preimages and `.flag` lists the (0-based) indices where Newton refinement failed to converge. `CrDiskMap::evalinv` and `AnnulusMap` return plain values instead, mirroring their MATLAB sources.
- **`StripMap`'s `ends` and `RectMap`'s `corners` are 1-indexed**, matching the MATLAB argument they translate — e.g. `StripMap(poly, {1, 4})`. This convention is deliberate throughout the port wherever an index is passed straight through from a `.m` source; see `CPP_PLAN.md` §3.5.
- **`AnnulusMap` takes scalars**, one point at a time, rather than vectors.
- **`Polygon` normalizes orientation.** Vertices are stored counterclockwise, reversing the input order if needed, exactly as `@polygon/polygon.m` does.

## Python bindings

NumPy-native bindings over the C++ port, built with [nanobind](https://github.com/wjakob/nanobind), live in [`python/`](python). Every map class, `AnnulusMap`, and `Polygon` are exposed with plain NumPy complex arrays in and out:

```python
import numpy as np
from sctoolbox import Polygon, DiskMap

verts = np.array([1j, -1 + 1j, -1 - 1j, 1 - 1j, 1, 0], dtype=complex)  # L-shape
m = DiskMap(Polygon(verts))

z = 0.6 * np.exp(1j * np.linspace(0, 2 * np.pi, 200))
w = m.eval(z)                       # disk -> polygon
print(np.max(np.abs(m.evalinv(w) - z)))   # round trip, ~1e-14
```

```sh
cd python && pip install .          # Eigen/nanobind fetched automatically
```

See [`python/README.md`](python/README.md) for the full API, the pytest suite, and a rendered gallery of conformal "carpet" plots for every map type.

## Tests and golden values

The C++ suite is a function-level comparison against MATLAB. Each `.m` function has recorded inputs and outputs in `tests/cpp/goldens_text/*.gold` (a plain-text format needing no MATLAB-file library), produced by generator scripts in `tests/cpp/generators/`. Regenerating them is only necessary if a `.m` source or a generator changes, and it does require MATLAB:

```matlab
addpath(genpath('/path/to/sc-toolbox'))
generateGoldens()      % -> tests/cpp/goldens/*.mat
exportGoldensText()    % -> tests/cpp/goldens_text/*.gold
```

The MATLAB toolbox has its own regression suite, independent of the port:

```matlab
result = sctool.runTests()
```

## Limitations

A few branches of the original are deliberately not ported and throw if reached — `crsplit`'s narrow-channel re-triangulation, `wquad1`'s line-segment continuation (whose `.m` source contains an indexing typo that would throw in MATLAB too), and `annulusmap`'s `'truncate'` option for unbounded outer polygons. Map constructors are restricted to the "polygon in, solve from an automatic initial guess" form. The `@crrectmap`, `@riesurfmap`, `@dscpolygons`, and `@scmapdiff` classes are not ported, and neither are any of the plotting or GUI routines.

`CPP_PLAN.md` documents the port in full: a phase-by-phase status table in §0, the design decisions per phase, and the numerical-parity notes — including several genuine bugs found in the 1998 MATLAB sources along the way.

## License

See [LICENSE](LICENSE). The original MATLAB toolbox is copyright Toby Driscoll.
