# C++ Translation Plan

Goal: exhaustive function-to-function, class-to-class translation of the SC Toolbox into C++, driven by fine-grained golden-value tests generated from the MATLAB reference implementation.

---

## 0. Current Status (read this first)

**All 13 phases are complete and passing.** As of the last run: **2242 assertions, 79 test cases, 0 failures.**

### Picking up on a new machine

Prerequisites: a C++17 compiler and CMake 3.16+. Eigen (3.4.0) and Catch2
(3.5.2) are fetched by `CMakeLists.txt` via `FetchContent`, pinned by version
and SHA256, so no system/Homebrew installs are needed and the build is
reproducible on any machine.

Build and test:

```sh
cmake -S cpp -B cpp/build
cmake --build cpp/build
./cpp/build/sctoolbox_tests           # or: ctest --test-dir cpp/build
```

`cpp/build/` is gitignored — it is regenerated from scratch by the commands above. Nothing else is needed to run the C++ tests: the goldens are checked in as plain text under `tests/cpp/goldens_text/*.gold`, and `CMakeLists.txt` points the test binary at that directory via the `GOLDENS_TEXT_DIR` compile definition. **MATLAB is not required to build or run the C++ tests** — only to regenerate goldens.

Regenerating goldens (MATLAB only, and only when a `.m` source or a generator changes):

```matlab
addpath(genpath('/path/to/sc-toolbox'))
generateGoldens()      % writes tests/cpp/goldens/*.mat  (calls rng(0) per group)
exportGoldensText()    % converts those to tests/cpp/goldens_text/*.gold
```

Both `goldens/*.mat` and `goldens_text/*.gold` are committed; the `.gold` text files are what the C++ side actually reads (§3.4).

### What is done

| Phase | Units | State |
|---|---|---|
| 1 | `gaussj`, `scqdata`, `scangle` | ✅ |
| 2 | `nesolve` family (`neqrdcmp`, `nechdcmp`, `nefdjac`, + internals) | ✅ |
| 3 | `scfix`, `scaddvtx`, `polygon`, `isinpoly` | ✅ |
| 4 | All `XXquad` (+ `stderiv`, `stquadh` pulled forward) | ✅ |
| 5 | All `XXderiv` | ✅ |
| 6 | All `XXmap` | ✅ |
| 7 | ODE RHS + custom Dormand-Prince `ode45`, `findz0` | ✅ |
| 8 | All `XXinvmap` | ✅ |
| 9 | `ellipkkp`, `ellipjc`, `r2strip`, `rcorners` | ✅ |
| 10 | All `XXparam` (incl. `crparam`'s triangulation pipeline) | ✅ |
| 11 | Map classes `DiskMap`/`HplMap`/`ExterMap`/`StripMap`/`RectMap`/`CrDiskMap` | ✅ |
| 12 | `AnnulusMap` (DSCPACK) | ✅ |
| 13 | `Moebius`, `Composite` | ✅ |

74 source files in `cpp/src/`, 51 test files in `cpp/tests/`, 21 golden groups in `tests/cpp/goldens_text/`.

### What is next

Nothing is outstanding in the translation plan. The port covers every function and class it set out to cover; what is left is either explicitly out of scope or a deliberately unported branch, both listed below.

### Explicitly out of scope

Never planned for translation, and no goldens exist: `@crrectmap`, `@riesurfmap`, `@dscpolygons`, `@scmapdiff`.

`@scmapinv` has no C++ class of its own, but its only purpose — an object whose `eval` is some map's `evalinv` — is what `Composite::inverseMember` produces, which is the only context `@scmapinv` is used in (see §12).

### Bug fixes found after the initial port

- `crtriang` edge ordering (`cpp/src/crtriang.cpp`): MATLAB stores each sub-polygon as a boolean mask and re-reads it with `find`, so its vertex list is always in increasing global order and `edge(:,enum)=idx(e)` is implicitly sorted. This port keeps a rotated vertex list (`i2` wraps around), so for polygons that subdivide more than once (n≥5) a stored diagonal edge's endpoints were left in position order, the base-case edge lookup missed it, and the resulting `-1` edge index caused an out-of-bounds `Eigen` write (a heap corruption that surfaced through the Python bindings). Fixed by sorting the stored endpoints by original vertex index; a defensive guard now throws instead of writing out of bounds. The MATLAB golden cases only cover quadrilaterals (one subdivision), so the path was never exercised — `test_crdiskmap.cpp` adds a golden-free self-consistency regression (octagon + L-shape round trip).

### Deliberately unported branches (throw `std::runtime_error` if reached)

- `crsplit`'s narrow-channel re-triangulation surgery — detection is ported, the mesh surgery is not (§9.6). Neither crdiskmap fixture triggers it.
- `wquad1`'s line-segment (`linearc==0`) continuation loop — the `.m` source has a self-indexing typo that would throw in MATLAB too (§11). Not exercised by any test case.
- `annulusmap`'s `'truncate'` option for unbounded outer polygons (`ishape==1`).
- Every map class's constructor is restricted to the "polygon in, parameter problem solved from an automatic initial guess" form; the MATLAB continuation-map / given-prevertex / tol-override constructor branches are not ported (§10).

Also skipped, but harmless: the pretty-printing methods of every class (`char`/`disp`/`display`), which have no numerical content.

### MATLAB-side changes this port depends on

These are committed alongside the C++ code and are required for golden generation:

- A `private_` static forwarding method on each of the 7 map classes (§3.2) — the only way to reach `private/` functions from a generator script.
- `@annulusmap/private/golden_fields.m` (§3.2) — reads internal fields past `annulusmap`'s restrictive `subsref`.
- `+sctool/nesolve.m`: `btrack = []` hoisted out of its `if (details(1) == 2)` guard, so the variable always exists before use.

---

## 1. MATLAB Built-ins That Need C++ Replacements

| MATLAB call | Where used | C++ replacement |
|---|---|---|
| `ode23(f, [0,0.5,1], z0, opts)` | `dinvmap`, `hpinvmap`, `deinvmap`, `stinvmap`, `rinvmap`, `crimap0` | Implement Bogacki-Shampine RK2(3) with dense output at requested t-points; or use Boost.Odeint `runge_kutta_fehlberg78` |
| `ode113(f, [0,0.5,1], z0, opts)` | `dinvmap` only | Substitute with `dopri5` (RK4(5)) from Boost.Odeint — the polish step corrects any extra error |
| `eig(tridiag)` | `sctool/gaussj.m` | LAPACK `dstev` (symmetric tridiagonal EVD) or `Eigen::SelfAdjointEigenSolver` |
| `gamma(x)` | `sctool/gaussj.m` | `std::tgamma(x)` (C++11) |
| `fsolve` | `@riesurfmap/private/rsparam.m` only | Port `nesolve` to C++ (used everywhere else); wire `rsparam` to the same C++ solver |
| `sort(z)` | scattered | `std::sort` |
| `norm(x, inf)` | `nesolve` family | `(x.cwiseAbs()).maxCoeff()` (Eigen) |
| Complex `abs`, `angle`, `real`, `imag`, `conj`, `exp`, `log`, `sqrt` | throughout | `std::abs`, `std::arg`, `std::real`, `std::imag`, `std::conj`, `std::exp`, `std::log`, `std::sqrt` on `std::complex<double>` |
| Matrix `\` (mldivide) | `nesolve` support routines have their own QR/Chol; isolated `\` elsewhere | `Eigen::ColPivHouseholderQR::solve` |
| `odeset('abstol', tol)` | inverse map calls | Pass tol as plain parameter to the ODE driver |

**No external MATLAB toolboxes are required** (fsolve/optimset in annulusmap and riesurfmap are either commented out or the only remaining use; nesolve covers all other nonlinear solving).

---

## 2. Recommended C++ Dependencies

- **Eigen 3** — dense linear algebra, complex vectors, matrix operations
- **Boost.Odeint** — ODE integration (substitute for ode23/ode113)
- **Catch2** (or GoogleTest) — unit test framework
- **LAPACK** (via Eigen's interface or directly) — symmetric tridiagonal eigenvalue in `gaussj`

---

## 3. Function-Level Golden-Value Test Infrastructure

### 3.1 Directory layout

```
tests/
  cpp/
    goldens/                   ← .mat files, one per functional unit
      gaussj.mat
      scqdata.mat
      scangle.mat
      scfix_disk.mat
      scfix_hpl.mat
      scfix_de.mat
      scfix_strip.mat
      scfix_rect.mat
      dquad.mat
      hpquad.mat
      dequad.mat
      stquad.mat
      crquad.mat
      dderiv.mat
      hpderiv.mat
      dederiv.mat
      stderiv.mat
      rderiv.mat
      crderiv.mat
      dmap.mat
      hpmap.mat
      demap.mat
      stmap.mat
      rmap.mat
      crmap.mat
      dimapfun.mat
      hpimapfun.mat
      stimapfun.mat
      rimapfun.mat
      dinvmap.mat
      hpinvmap.mat
      deinvmap.mat
      stinvmap.mat
      rinvmap.mat
      crinvmap.mat
      dparam.mat
      hpparam.mat
      deparam.mat
      stparam.mat
      rparam.mat
      crparam.mat
      ellipkkp.mat
      ellipjc.mat
      r2strip.mat
      neqrdcmp.mat
      nechdcmp.mat
      nefdjac.mat
      nemodel.mat
      nelnsrch.mat
      nehook.mat
      nesolve.mat
      annulus_qinit.mat
      annulus_wquad.mat
      annulus_dscsolv.mat
      annulus_eval.mat
    generateGoldens.m          ← master script; calls all gen_*.m helpers
    generators/
      gen_gaussj.m
      gen_scqdata.m
      gen_scangle_scfix.m
      gen_diskmap_private.m    ← dquad, dderiv, dmap, dimapfun, dinvmap, dparam
      gen_hplmap_private.m
      gen_extermap_private.m
      gen_stripmap_private.m
      gen_rectmap_private.m    ← includes ellipkkp, ellipjc, r2strip
      gen_crdiskmap_private.m
      gen_annulus_private.m
      gen_nesolve.m
    testCppGoldens.m           ← runs verification: recomputes and diffs vs goldens
```

### 3.2 How to access private functions from generator scripts

MATLAB private directories are only accessible from `.m` files that live directly inside the owning class folder (including its own `private/` subfolder) — `addpath` cannot be used on any folder named `private` (MATLAB blocks this outright), so the directly-add-to-path approach does not work.

The mechanism that does work: each map class (`diskmap`, `hplmap`, `extermap`, `stripmap`, `rectmap`, `crdiskmap`, `annulusmap`) has a generic Static forwarding method added to its classdef file:

```matlab
methods (Static)
    function varargout = private_(name, varargin)
        [varargout{1:nargout}] = feval(name, varargin{:});
    end
end
```

Because `private_` is itself defined inside the class folder, MATLAB resolves the `feval(name, ...)` call against that folder's `private/` functions. Generator scripts then call private functions like:

```matlab
[z, c, qdat] = diskmap.private_('dparam', w, beta, [], opt);
I            = diskmap.private_('dquad', z1, z2, sing1, z, beta, qdat);
```

This is the only pattern used throughout `tests/cpp/generators/`. When writing the corresponding signature, always read the actual `function ... = name(...)` line in `@ClassName/private/name.m` rather than assuming an argument order — several generators in this codebase were originally written with wrong argument orders (e.g. `hpmap`/`demap`/`rmap`/`stmap` all take `w` before `z`; `rinvmap` takes `L` before `qdat`) that were only caught by cross-checking the real signatures.

**`annulusmap` is a special case.** Unlike the `scmap`-derived classes (which inherit `@scmap/subsref.m`'s safe fallback — only `()`-indexing is special-cased, so ordinary dot-notation property access like `m.qdata` works fine from anywhere), `annulusmap` has its own restrictive `subsref.m` that throws an error for ANY dot-notation access (`m.M`, `m.u`, etc.) performed by code outside the class folder — and this restriction applies even to `builtin('subsref', ...)` calls and even when invoked indirectly via `feval` through `private_` from outside, because MATLAB's virtual dispatch for the overridden `subsref` happens regardless of call-site indirection. The only way to read internal fields of an `annulusmap` object from a generator script is to add a small accessor function *inside* the class's own folder — code that "ships with the class" gets full internal access regardless of how it's invoked, but code outside the folder never does. This codebase includes such a helper, `@annulusmap/private/golden_fields.m`, which bundles the needed fields into a returned struct:

```matlab
f = annulusmap.private_('golden_fields', m);
% f.M, f.N, f.ALFA0, f.ALFA1, f.u, f.c, f.w0, f.w1, f.phi0, f.phi1
```

Use this `golden_fields.m` pattern as the template for any other class encountered later that turns out to have a similarly restrictive `subsref`.

Some classes' `get(m, 'propname')` accessor methods (e.g. `@extermap/get.m`, `@rectmap/get.m`) only recognize a small set of property-name prefixes and silently return `[]` with a `warning(...)` (not an error) for anything else — e.g. `get(m,'qdata')` on an `extermap`, or `get(m,'rectangle')` on a `rectmap`. This can silently corrupt downstream golden values without any visible failure. Prefer direct dot-notation property access (`m.qdata`, `m.stripL`, `m.prevertex`, `m.constant`) over `get(...)` for any class other than `annulusmap`.

### 3.3 Golden file schema

Each `.mat` file stores one struct array `cases`, where every element has:

```
cases(k).desc     — string describing the test case
cases(k).inputs   — struct of named inputs (matching the MATLAB function signature)
cases(k).outputs  — struct of named outputs
cases(k).tol      — absolute tolerance for C++ comparison (typically 1e-11)
```

Use multiple cases per file to cover: interior points, boundary-adjacent points, singular prevertices, near-degenerate polygons, polygons with infinite vertices.

### 3.4 `.gold` text export (actual interchange format used by the C++ side)

JSON was considered and rejected: `jsonencode` has no native complex-number or
Inf/NaN support without hand-rolled conversion, and reading `.mat` directly
from C++ would require `libmatio`. Instead `tests/cpp/exportGoldensText.m`
converts any `tests/cpp/goldens/<name>.mat` into a plain-text
`tests/cpp/goldens_text/<name>.gold`:

```
GROUP <groupname>
CASE
DESC <description>
TOL <tolerance>
INPUT <field> <rows> <cols> R|C
<row-major values, one row per line, "%.17g" precision>
...
OUTPUT <field> <rows> <cols> R|C
...
ENDCASE
...
ENDGROUP
```

`R` fields have one value per element per line; `C` fields have two
(`real imag`). String fields use `INPUT/OUTPUT <field> STR` followed by a
literal line. Non-numeric, non-char fields (e.g. `crdiskmap`'s `Q` qlgraph
object) are written as `INPUT/OUTPUT <field> SKIP` and the exporter emits a
MATLAB `warning(...)` — these tiers will need a per-field strategy when
reached (not yet designed). `Inf`/`-Inf`/`NaN` round-trip natively since
`std::stod`/`strtod` parse them per C99/C++11. The C++ side parses this with
the self-contained reader in `cpp/tests/golden_reader.hpp`
(`golden::loadGoldens` → `GroupMap` of `golden::Case`, each with `.in(name)` /
`.out(name)` accessors returning `Eigen::MatrixXcd`).

### 3.5 Actual C++ project layout (Phase 1 / Tier 0 and Phase 2 complete)

```
cpp/
  CMakeLists.txt              ← Eigen3 + Catch2 via find_package, ctest integration
  include/sctoolbox/
    gaussj.hpp
    scqdata.hpp
    scangle.hpp
    neqrdcmp.hpp
    nechdcmp.hpp
    nefdjac.hpp
    nesolve.hpp
  src/
    gaussj.cpp
    scqdata.cpp
    scangle.cpp
    neqrdcmp.cpp
    nechdcmp.cpp
    nefdjac.cpp
    nesolve.cpp
  tests/
    golden_reader.hpp
    test_gaussj.cpp
    test_scqdata.cpp
    test_scangle.cpp
    test_neqrdcmp.cpp
    test_nechdcmp.cpp
    test_nefdjac.cpp
    test_nesolve.cpp
```

Build: `cmake -S cpp -B cpp/build -DCMAKE_PREFIX_PATH=/opt/homebrew && cmake --build cpp/build`, run with `cpp/build/sctoolbox_tests`. As of this writing all functions across Phase 1-4 and the unblocked portions of Phase 5-6 pass against MATLAB-generated goldens (1124 assertions, 23 test cases, 0 failures). `nemodel`, `nelnsrch`, `nehook`, `netrust`, `nebroyuf`, `neinck`, `nersolv`, `neqrsolv`, `neconest`, and `nefn` are ported as anonymous-namespace internals of `nesolve.cpp` (not separately exposed/tested) and are exercised indirectly through `nesolve`'s residual-based golden comparisons rather than function-level goldens of their own.

`sctoolbox::isinpoly` (`cpp/src/isinpoly.cpp`) ports `+sctool/isinpoly.m`'s argument-principle winding test, including its complex `sign()` (`csign`, unit-modulus normalization, 0 at 0) and its two distinct round-then-divide vs. divide-then-round formulas for the no-boundary-hit vs. boundary-hit cases. `sctoolbox::Polygon` (`cpp/src/polygon.cpp`) ports `@polygon/polygon.m` + `@polygon/angle.m` as a single constructor path: angle computation (with the all-collinear early-return special case preserved), the orientation check/reversal, and the multiple-sheeted warning (not fatal, so omitted as a no-op). Both are covered by `test_isinpoly.cpp`/`isinpoly.gold` and `test_polygon.cpp`/`polygon.gold`, generated via `gen_isinpoly.m`/`gen_polygon.m`.

Phase 4 (`dquad.cpp`, `hpquad.cpp`, `dequad.cpp`, `stquad.cpp`, `crquad.cpp`) all share the same `sing1`-as-1-indexed-MATLAB-convention interface established for `scfix`: a `0` sentinel means "no singularity at this endpoint," otherwise the value is the 1-indexed position within `z`, converted to a 0-indexed row only at the point of array access. Each function mirrors its `.m` source almost line-for-line, including the qdat column-selection formula `ind = sng if sng>=1 else (column n)` (originally `rem(sng+n,n+1)+1` in MATLAB) and the adaptive-subdivision `while (dist < 1)` loop. `stquad.cpp` depends on a newly-ported `stderiv.cpp` (the Tier-2 derivative function `@stripmap/private/stderiv.m`, pulled forward because `stquad.m` calls it directly as its integrand) — `stderiv` strips the strip map's two infinite "end" prevertices out of `z`/`beta` internally and re-indexes the optional Gauss-Jacobi normalization index `j` around the removal, exactly as MATLAB does. `crquad.cpp` differs from the other four in structure (integrates from `z1(k)` to an implicit right endpoint of 0, via geometric dyadic-panel subdivision `panels = max(1, ceil(-log(mindist)/log(2)))` rather than the single-shrink `while` loop) and additionally masks out vertices with `beta≈0` (`ignore`) from the product, re-indexing the singularity position into the masked `keep` list — ported faithfully including the quirk that a detected coincident-prevertex panel resets the entire accumulated integral `I(k)` to 0 (not just that panel's contribution).

Phase 5 ports the remaining self-contained `XXderiv` functions: `dderiv.cpp` (`f'(zp) = c * prod_i (1 - zp/z(i))^beta(i)`), `hpderiv.cpp` (same log-sum-product form but with `terms = zp - z(i)` over only the finite prevertices, infinite entries of `z`/`beta` stripped first), and `dederiv.cpp` (disk-product form with an extra appended `-2*log(zp)` term accounting for the exterior map's implicit singularity at the origin). `stderiv` was already ported in Phase 4 (pulled forward as a `stquad` dependency). Two derivative functions remain genuinely blocked and are not yet portable: `rderiv` (`@rectmap/private/rderiv.m`) depends on `r2strip`/`ellipkkp`/`ellipjc` (Phase 9, Landen-transformation elliptic functions) and `stripmap.deriv`; `crderiv` (`@crdiskmap/private/crderiv.m`) depends on `crembed`/`crspread` (cross-ratio quadrilateral-graph embedding machinery built during `crparam`, Phase 10). Both will be revisited once their respective phases land.

Phase 6 ports the four self-contained `XXmap` forward-evaluation functions: `dmap.cpp`, `hpmap.cpp`, `demap.cpp`, `stmap.cpp`. Each finds, for every target point, the nearest prevertex (screening out exact vertex images and, for `hpmap`/`stmap`, points already at infinity), then integrates the relevant `XXquad`/`XXquadh` from that prevertex to the target, handling the "bad point" case (closest prevertex maps to an infinite polygon vertex, so integration must instead start from/route through a neighboring finite prevertex or a conformal-center basis point) per-map-type exactly as the `.m` source does. `stmap.cpp` additionally required porting `stquadh.cpp` (`@stripmap/private/stquadh.m`), the adaptive horizontal-interval subdivider used for the "bad point" routing segment, which recursively splits the integration interval until no singularity violates the alpha-rule horizontal/vertical distance bound, then delegates to `stquad` per safe sub-interval. While porting `dmap`, a latent bug was found and fixed in `tests/cpp/generators/gen_diskmap_private.m`: the `cases_dmap`/`cases_dinvmap` blocks called `dmap` with the `w` and `z` arguments swapped and `z` passed as an empty array, which degenerately produced all-zero golden outputs (a weak, non-discriminating test) rather than erroring — `diskmap_private.gold` was regenerated after the fix to produce meaningful non-trivial golden values. `rmap` (`@rectmap/private/rmap.m`) and `crmap` (`@crdiskmap/private/crmap.m`) remain blocked on the same dependencies as `rderiv`/`crderiv` (`r2strip`/elliptic functions, and `crembed`/`crspread`, respectively).

`scfix.cpp`/`scaddvtx.cpp` mirror the MATLAB source's 1-indexed arithmetic directly (permutation helpers take/return 1-indexed position lists) rather than rewriting everything to 0-indexed form up front — this was deliberate, since scfix.m's renumber/insert logic is dense enough that re-deriving the index algebra from scratch in 0-indexed terms is a major source of off-by-one risk; only the final array element accesses drop to 0-indexed. `golden_reader.hpp`'s `STR` fields (e.g. scfix's `type` argument) are now captured into `Case::inputStrings`/`outputStrings` rather than discarded, since `test_scfix.cpp` needs the type string per case.

`golden_reader.hpp`'s `readMatrix` parses each numeric token with `strtod` rather than `std::istream::operator>>(double&)`, since the latter fails to parse "NaN"/"Inf" tokens under libc++ on this platform — this matters because several Phase 2 goldens (`nefdjac`'s degenerate zero-step-size case) legitimately contain non-finite values.

---

## 4. Functional Units and Test Cases — Tier by Tier

### Tier 0 — Pure Math Primitives

#### `gaussj(n, alf, bet)` → `[z, w]`

Cases:
- `(4, 0, 0)` — standard Gauss-Legendre (known nodes and weights)
- `(8, 0, 0)` — 8-point Gauss-Legendre
- `(8, 0, 0.5)` — non-integer beta typical of SC
- `(8, 0, -0.5)` — negative beta (π-angle vertex)
- `(8, 0, 1.5)` — large beta (reentrant corner)
- `(1, 0, 0)` — edge case n=1

Verify: `sum(w) == 2^(alf+bet+1)*B(alf+1, bet+1)`, nodes in `(-1, 1)`, monotonically increasing.

#### `scqdata(beta, nqpts)` → `qdat`

Cases:
- Square polygon `beta = [0.5, 0.5, 0.5, 0.5]`, `nqpts = 8`
- L-shaped polygon with a reentrant corner `beta = [0.5,0.5,1.5,0.5,0.5,0.5]`, `nqpts = 8`
- Polygon with an infinite vertex (`beta(j) < -1` skipped in `scqdata`)
- `nqpts = 4` and `nqpts = 12`

Verify: matrix dimensions `(nqpts, 2*(n+1))`, column pairs sum correctly for each Jacobi weight.

#### `sctool.scangle(w)` → `beta`

Cases:
- Unit square vertices (expect all 0.5)
- Triangle (expect all 1/3)
- L-shape
- Pentagon
- Polygon with an infinite vertex

---

### Tier 1 — Quadrature

All quadrature functions share the same interface pattern: `XXquad(z1, z2, sing1, z, beta, qdat)`.  
For each, generate cases by:
1. Fixing a small solved map (pentagon, square, L-shape from `fixturePolygons`)
2. Selecting representative integration segments: singularity-to-interior, interior-to-interior, singularity-to-singularity (requires split)
3. Comparing against a finer-nqpts reference evaluation

#### `dquad(z1, z2, sing1, z, beta, qdat)` — disk

Cases:
- `sing1 = 1` (integrate from prevertex 1 to mid-arc)
- `sing1 = 0` (regular interval)
- Multiple simultaneous intervals (vectorized)
- Interval requiring adaptive subdivision (singularity at distance < 0.5 of interval length)

#### `hpquad(z1, z2, sing1, z, beta, qdat)` — half-plane

Same structure as `dquad` but prevertices on real line; include `Inf` prevertex cases.

#### `dequad(z1, z2, sing1, z, beta, qdat)` — exterior disk

Include the `sing1 = 0` case involving the origin (extra singularity for exterior maps).

#### `stquad(z1, z2, sing1, z, beta, qdat, L)` — strip

Include the exponential coordinate transformation that maps strip to half-plane.

#### `crquad(...)` — cross-ratio disk

Use the quadrilateral test polygons from `fixturePolygons`.

---

### Tier 2 — Derivatives

Pattern: `XXderiv(zp, z, beta, c)` computes `f'(zp) = c * ∏(zp - z_k)^beta_k`.

#### `dderiv(zp, z, beta, c)` — disk
- Points well inside disk
- Points near (but not at) prevertices
- Multiple simultaneous points
- Verify against finite difference of `dmap`: `|dderiv(z) - (dmap(z+h) - dmap(z-h))/(2h)| < 1e-5` for `h=1e-6`

#### `hpderiv`, `dederiv`, `stderiv`, `rderiv`, `crderiv`

Same structure. For `rderiv`, the derivative involves elliptic functions (via `r2strip`).

---

### Tier 3 — Forward Maps

`XXmap(zp, z, c, L, beta, qdat)` integrates from a reference prevertex to `zp`.

#### `dmap(zp, z, c, L, beta, qdat)` — disk
- Single point
- Column vector of points
- Point at prevertex (returns corresponding polygon vertex)
- Point exactly at `z(1)` (boundary singularity)
- Verify: `dmap(z(k)) ≈ w(k)` for all prevertices k

#### `hpmap`, `demap`, `stmap`, `rmap`, `crmap`

Same verification strategy. For `rmap`, additionally verify that the four corners map to the four polygon vertices.

---

### Tier 4 — ODE Right-Hand Sides

`XXimapfun(t, y, fdat)` is the ODE integrand for inverse mapping.  
The argument `y` packs `[real(z); imag(z)]` into a real vector.

Cases for each:
- Single trajectory point at `t = 0, 0.5, 1`
- Multiple packed trajectories
- Verify: RHS is `1/f'(z(t))` scaled by the chord direction `(w_target - w_0)`

These are tested by integrating a known path and comparing the ODE solution to the known preimage.

---

### Tier 5 — Inverse Maps

`XXinvmap(wp, w, beta, z, c, qdat)` returns preimage points.

Cases (use solved fixture maps):
- Single interior point
- Grid of interior points
- Point at a polygon vertex (should return corresponding prevertex)
- Point near a corner (numerically challenging)
- Batch of points (vectorized)

Round-trip check: `eval(m, evalinv(m, w)) ≈ w` to `1e-10`.

---

### Tier 6 — Nonlinear Solver Components

Test each sub-routine of `nesolve` independently.

#### `neqrdcmp(A)` — Householder QR
- Tall, square, and nearly singular matrices
- Verify: `Q * R ≈ A` where `Q` is reconstructed from Householder reflectors
- Verify: `R` upper triangular, diagonal non-negative

#### `nechdcmp(H, macheps)` — Perturbed Cholesky
- Positive definite input (perturbation should be zero)
- Indefinite input (diagonal augmentation kicks in)
- Verify: `L * L' ≈ H + mu*I`

#### `nefdjac(fvec, Fvec, x, Sx, details, nofun, fparam)` — Finite-diff Jacobian
- Small test function `f(x) = [x(1)^2 + x(2); x(1) - x(2)^2]`
- Compare to analytic Jacobian to `1e-6`

#### `nemodel(F, J, g, SF, Sx, iflg)` — Linear model
- Verify `sN = -J\F` (Newton step) for well-conditioned system

#### `nelnsrch` — Line search
- Starting from non-solution, verify Armijo condition is satisfied at returned step

#### `nehook` — Trust-region step
- Verify step length ≤ trust radius
- Verify step solves `(J'J + mu*I) s = -J'F` for some `mu ≥ 0`

#### `nesolve` — Full solver
- Simple 2D test system `f(x) = [x(1)^2 - 1; x(2)^2 - 4]`, start `x0 = [2; 3]`
- Ill-conditioned system to test Broyden fallback
- System that requires trust-region globalization
- Verify `termcode == 1` (normal termination) and `norm(F(xf), inf) < tol`

---

### Tier 7 — Parameter Problems

These are the core scientific computations. Inputs are polygon vertices `w` and turning angles `beta`; outputs are prevertices `z`, constant `c`, and quadrature data.

#### `dparam(w, beta, z0, opt)` — disk
Cases (use polygons from `fixturePolygons`):
- Square
- Pentagon
- L-shape (reentrant corner)
- Polygon with one infinite side
- Provide initial guess `z0` (warm start)

Verify:
- `sum(beta) == -2` (closed polygon constraint)
- `norm(dpfun(z, w, beta, qdat), inf) < tol` (residual of parameter equations)
- Round-trip: `dmap(z, z, c, [], beta, qdat) ≈ w`

#### `hpparam`, `deparam`, `stparam`, `rparam`, `crparam`

Same verification strategy with domain-appropriate polygons and residual functions.

---

### Tier 8 — Elliptic Functions (rectmap)

#### `ellipkkp(L)` → `[K, Kp]`
- Known values: `L = 1/sqrt(2)` → `K = Kp = Γ(1/4)^2 / (4*sqrt(π))` ≈ 1.8541
- Range of `L` in `(0, 1)`: 0.1, 0.3, 0.5, 0.7, 0.9, 0.99
- Verify: `K * Kp` relation, AGM convergence to `< 1e-14`

#### `ellipjc(u, L)` → `[sn, cn, dn]`
- Real argument: compare to MATLAB's `ellipj`
- Pure imaginary argument (exercises the complex addition formula)
- General complex argument
- Verify: `sn^2 + cn^2 == 1`, `dn^2 + L^2 * sn^2 == 1`

#### `r2strip(zp, z, c, L)` — rectangle to strip conversion
- Verify it matches `stmap(ellipjc(...))` composition
- Round-trip through `r2strip` and its inverse

---

### Tier 9 — Annulus Private Functions

#### `qinit(map, nptq)` → `qwork`
- Verify dimensions match `(M + N)` vertices and `nptq` quadrature points

#### `wquad(...)` and `wquad1(...)`
- Fix a solved annulus map; verify quadrature values match finer reference

#### `dscsolv(iguess, nptq, qwork, isUnbounded, linearc, map)` → `[u, c, w0, w1, phi0, phi1]`
- Use square-annulus test case from `testAnnulus.m`
- Verify: round-trip `annulusmap eval → evalinv` to `1e-10`

#### `annulusmap eval / evalinv`
- Interior points on annular grid
- Points on outer and inner boundary
- Round-trip check

---

### Tier 10 — Polygon Class

#### Constructor and accessors
```
polygon([0,1,1+1i,1i])                   % unit square
polygon([0,1,1+1i,1i], [.5,.5,.5,.5])   % with explicit angles
vertex(p), angle(p), length(p), isinf(p)
```

#### `scfix`
- Verify cyclic renumbering for disk constraints
- Verify vertex insertion when needed
- All 5 types: `'d'`, `'hp'`, `'de'`, `'st'`, `'r'`

---

## 5. Golden Generation Script Skeleton

```matlab
% tests/cpp/generateGoldens.m
function generateGoldens(which)
ALL = {'gaussj','scqdata','scangle_scfix', ...
       'diskmap_private','hplmap_private','extermap_private', ...
       'stripmap_private','rectmap_private','crdiskmap_private', ...
       'annulus_private','nesolve'};
if nargin < 1, which = ALL; end
outDir = fullfile(fileparts(mfilename('fullpath')), 'goldens');
if ~exist(outDir,'dir'), mkdir(outDir); end
for i = 1:numel(which)
    genFn = str2func(['gen_' which{i}]);
    fprintf('--- %s ---\n', which{i});
    genFn(outDir);
end
fprintf('Done.\n');
end
```

Each `gen_XXX.m` function follows the pattern:

```matlab
% tests/cpp/generators/gen_gaussj.m
function gen_gaussj(outDir)
cases = struct('desc',{},'inputs',{},'outputs',{},'tol',{});

% Case 1: Gauss-Legendre (alf=bet=0)
[z,w] = sctool.gaussj(8, 0, 0);
cases(end+1) = struct( ...
    'desc', 'n=8 alf=0 bet=0 (Gauss-Legendre)', ...
    'inputs',  struct('n',8,'alf',0,'bet',0), ...
    'outputs', struct('z',z,'w',w), ...
    'tol', 1e-14);

% Case 2: Jacobi with typical SC exponent
[z,w] = sctool.gaussj(8, 0, 0.5);
cases(end+1) = struct( ...
    'desc', 'n=8 alf=0 bet=0.5', ...
    'inputs',  struct('n',8,'alf',0,'bet',0.5), ...
    'outputs', struct('z',z,'w',w), ...
    'tol', 1e-14);

% ... more cases ...

save(fullfile(outDir,'gaussj.mat'), 'cases');
fprintf('  gaussj: %d cases\n', numel(cases));
end
```

---

## 6. C++ Test Verification Script

```matlab
% tests/cpp/testCppGoldens.m
% Reload all goldens and recompute to verify nothing has drifted.
% Also serves as the acceptance spec for the C++ port.
function testCppGoldens()
goldenDir = fullfile(fileparts(mfilename('fullpath')), 'goldens');
files = dir(fullfile(goldenDir, '*.mat'));
nPass = 0; nFail = 0;
for i = 1:numel(files)
    s = load(fullfile(goldenDir, files(i).name));
    % dispatch to per-type verifier ...
end
fprintf('%d passed, %d failed\n', nPass, nFail);
end
```

---

## 7. C++ Translation Order

Translate in dependency order. Each tier can only begin after the previous tier's C++ functions pass their golden tests.

| Phase | Functional units | Key C++ decisions |
|---|---|---|
| 1 ✅ | `gaussj`, `scqdata`, `scangle` | Use `Eigen::SelfAdjointEigenSolver` for tridiagonal EVD; `std::tgamma` |
| 2 ✅ | `neqrdcmp`, `nechdcmp`, `nefdjac`, `nemodel`, `nelnsrch`, `nehook`, `nesolve` | Port Dennis & Schnabel verbatim into C++; keep the same details-vector interface or refactor to a struct |
| 3 ✅ | `scfix`, `polygon` class | C++ `Polygon` class with `vertices` and `angles` as `Eigen::VectorXcd` / `Eigen::VectorXd`; `scfix` as free function |
| 4 ✅ | All `XXquad` functions | One quadrature template; specialize for each domain |
| 5 ✅ | `dderiv`, `hpderiv`, `dederiv`, `stderiv`, `rderiv`, `crderiv` | Thin wrappers around the log-sum product formula; `crderiv` additionally needs `crembed`/`crspread`/`moebius` |
| 6 ✅ | `dmap`, `hpmap`, `demap`, `stmap`, `rmap`, `crmap` | Integrate `XXquad` from reference prevertex |
| 7 ✅ | `XXimapfun` (ODE RHS) | Custom Dormand-Prince `ode45` (no Boost dependency); one ODE-RHS lambda per map type, fed by `findz0` for initial guesses |
| 8 ✅ | `XXinvmap` (inverse maps) | `dinvmap`/`hpinvmap`/`deinvmap`/`stinvmap`/`rinvmap`/`crinvmap` — ODE continuation + Newton polish, matching each map's boundary reflection rule; `rinvmap` additionally clamps to the rectangle bounding box (`rectproject`) instead of reflecting; `crinvmap` uses `isinpoly` to pick a starting quadrilateral, then `crimap0` + `crgather` |
| 9 ✅ | `ellipkkp`, `ellipjc`, `r2strip`, `rcorners` | Port Landen-transformation algorithm from `ellipjc.m`; `findz0` extended with the `"r"` (rectangle) prefix |
| 10 ✅ | `dparam`, `hpparam`, `deparam`, `stparam`, `rparam`, `crparam` | Wire `nesolve` + domain residual functions; `dparam`/`deparam` additionally needed `sctool.dabsquad` (ported alongside `dparam`); `stparam` reuses `stquad`/`stquadh` and renumbers around the strip's two end-vertices; `rparam` additionally needed `rptrnsfm` (Trefethen-style short-edge transform) and a post-`nesolve` Newton refinement onto the rectangle boundary via `r2strip`; `crparam` needed a full polygon-triangulation pipeline (`crtriang`, `crcdt`, `crqgraph`, `crsplit`, `crossrat`, `craffine`, `crfixwc`) plus `nesolve`'s identity-initial-Jacobian variant (`nesolvei.m`) |
| 11 ✅ | Map classes (`diskmap`, `hplmap`, `extermap`, `stripmap`, `rectmap`, `crdiskmap`) | C++ classes with `eval`, `evalinv`, `evaldiff`, `accuracy`; constructors call the matching `XXparam` directly (not via a from-scratch port of every MATLAB constructor branch) |
| 12 ✅ | `annulusmap` | Standalone; depends on its own quadrature and `nesolve` |
| 13 ✅ | `moebius`, `composite` | `Moebius` ports the three-point constructor's infinity branches verbatim; `Composite` replaces MATLAB's duck typing with a forward/inverse callable pair per member, which is what lets `inverse()` reverse-and-swap the chain the way `@composite/inv.m` does |

---

## 8. Naming and File Layout (proposed C++ repo structure)

```
sc-toolbox-cpp/
  include/sctoolbox/
    polygon.hpp
    scmapopt.hpp
    quadrature.hpp       ← gaussj, scqdata
    scfix.hpp
    nesolve/             ← all nesolve support headers
    maps/
      scmap.hpp          ← base class
      diskmap.hpp
      hplmap.hpp
      extermap.hpp
      stripmap.hpp
      rectmap.hpp
      crdiskmap.hpp
      annulusmap.hpp
      moebius.hpp
      composite.hpp
  src/
    (*.cpp implementations)
  tests/
    (Catch2 test files, one per functional unit mirroring the MATLAB golden files)
    goldens/             ← symlink or copy of tests/cpp/goldens/ from MATLAB repo
```

---

## 9. Numerical Parity Requirements

The C++ port must match MATLAB to within `5e-11` (absolute) on all golden test cases. This matches the tolerance used by the existing MATLAB regression tests (`fixture.tol * 50` where `tol = 1e-12`).

The only known sources of deviation:
- ODE solver substitution (`ode113` → `dopri5`): Newton polish brings both to `< 1e-10`, so this is acceptable
- `gamma` vs `std::tgamma`: identical to IEEE double precision
- Tridiagonal EVD (`eig` vs `Eigen::SelfAdjointEigenSolver`): both are backward stable; results will be bit-for-bit compatible

### 9.1 `sctool.findz0` draws from the global RNG

`+sctool/findz0.m` (used by `deinvmap`, `crinvmap`, and other inverse-map private functions to pick an initial point/search direction for ODE continuation) calls `rand(1)` to decide whether to "abandon midpoints" during its search. This means deinvmap-family results are **not purely a function of their documented arguments** — they also depend on the ambient global RNG state, which in turn depends on what ran earlier in the same MATLAB session/process.

This was discovered when `gen_extermap_private.m`'s `deinvmap` golden cases produced different "Check solution; maximum residual = ..." warnings (and slightly different `zp` results) depending on whether the generator was run standalone vs. as part of the full `generateGoldens()` batch. `generateGoldens.m` now calls `rng(0)` before invoking each generator group, so results are reproducible across runs regardless of group ordering — but this also means:
- The recorded golden `zp_inv` values for `deinvmap`/`crinvmap`-family cases are tied to MATLAB's specific RNG algorithm and the exact sequence/count of `rand` calls made by `findz0` for that case. They are not a "pure" mathematical golden value.
- The C++ port does not need to replicate MATLAB's RNG bit-for-bit — `findz0`'s `rand(1)` branch is a heuristic fallback for ODE search robustness, not something that needs numerical parity. When porting `findz0`, a deterministic substitute (e.g. always taking the same branch, or a fixed pseudo-random sequence) is acceptable; just don't expect the C++ `deinvmap`/`crinvmap` outputs to match these specific MATLAB goldens bit-for-bit if `findz0`'s branch choice differs. Prefer comparing such cases by residual (`|wp - forwardmap(zp)| < tol`) rather than exact `zp` equality.

**Confirmed in practice (extermap `deinvmap`):** the original `gen_extermap_private.m` golden case fed `deinvmap` points whose forward image (`wp`) sits near a dense cluster of `beta=-0.5` prevertices on the unit circle. Verified directly against MATLAB (varying `tol`/`maxiter` via `scinvopt`-style options) that *MATLAB's own* `deinvmap` does not converge for these inputs — the residual stalls around 0.02-0.3 regardless of iteration budget. This is a basin/conditioning issue in the algorithm itself, not a `demap`/`dequad` bug. Fixed by changing that golden case's input points to well-conditioned ones close to (but outside) the unit circle, with a tightened `[0, 1e-12, 80]` options vector, so the recorded `zp` is an accurately-converged round-trip value rather than MATLAB's own non-convergent artifact. The C++ test now passes the same tightened options.

### 9.2 `rmap`/`rinvmap`'s augmented qdat needs an extra "neutral" column, not just the two spliced ones

`rmap.m`/`rderiv.m` insert two synthetic `Inf`/`-Inf` "prevertices" into the strip-mapped `zs`/`ws`/`bs` arrays (length goes from `n` to `n+2`) so they can be fed through `stripmap`'s `stmap`/`stderiv`, which require their prevertex list to literally contain the strip's infinite ends. The qdat columns must be re-spliced to match: MATLAB does this with `idx = [1:ends(1) n+1 ends(1)+1:ends(2) n+1 ends(2)+1:n n+1]` (length `n+3`, **three** copies of the filler column `n+1`, not two) and `qdat = qdat(:,[idx idx+n+1])`.

It is tempting (and was an initial mistake in the C++ port) to assume the third `n+1` filler is dead/unused, since only two `Inf`/`-Inf` slots were inserted into `zs`/`ws`/`bs`. It is not dead: `stmap`/`stquad`'s internal column indexing assumes qdat is laid out in the same "(prevertex count)+1"-sized-block convention that `scqdata` itself produces (one real node/weight column per prevertex, plus one trailing "neutral" `beta=0` column used whenever a quadrature segment has no local singularity). Since the augmented prevertex count is `n+2`, the augmented qdat needs `n+2+1 = n+3` column-pairs to supply that same trailing neutral column — which is exactly what MATLAB's third filler provides. Dropping it (producing only `2*(n+2)` columns instead of `2*(n+3)`) causes `stmap`/`stquad` to index one column short of the neutral column whenever a quadrature node has no singularity, silently reading the wrong (real-data) column and producing wrong — not NaN/crashing — results. Caught by reproducing MATLAB's exact intermediate `qdataug` matrix and comparing column-by-column; fixed in `rectAugQdat()` (`cpp/src/rectmap_internal.cpp`) by building the full `n+3`-length `idx` array (with all three fillers) rather than `n+2`.

### 9.3 `crdiskmap`'s quadrilateral graph `Q` is plain numeric data — it was just never exported

`crderiv`/`crmap`/`crinvmap` (`@crdiskmap/private/*.m`) take a `Q` argument (the "quadrilateral graph" built by `crqgraph.m` from a Delaunay-ish triangulation, deferred along with the rest of `crparam` to the not-yet-started Phase 10). §3.4 originally flagged `Q` as a non-numeric MATLAB struct that the `.gold` exporter would have to `SKIP` with a "needs a per-field strategy when reached" warning. In practice `Q` turned out to be three plain numeric/logical matrices (`Q.qlvert` 4×n3, `Q.qledge` 4×n3, `Q.adjacent` n3×n3) — `crqgraph.m`'s only non-trivial output, with no actual object/handle involved. `gen_crdiskmap_private.m` now flattens these into top-level `Qqlvert`/`Qqledge`/`Qadjacent` golden fields (1-indexed, like everything else MATLAB emits), which the existing numeric-only exporter already serializes correctly — no exporter change was needed, and `crqgraph`/`crtriang`/`crcdt` (the triangulation that builds `Q` from scratch) remain unported, exactly mirroring how `rderiv`/`rmap`/`rinvmap` take `z`/`c`/`L` directly instead of requiring `rparam`. The C++ side represents this as a plain `QGraph` struct (`cpp/include/sctoolbox/qgraph.hpp`) of 0-indexed `Eigen::MatrixXi`s.

**Latent MATLAB bug found in the process:** `crderiv.m`/`crmap.m` both compute `[~,idx] = min(abs(zl))` to pick, per query point, the best-conditioned embedding out of `n3` quadrilaterals. When a polygon has only one quadrilateral (`n3==1`, e.g. a plain unit square) and `zl` is the resulting `1×m` matrix, MATLAB's `min` treats a 1-row matrix as a *vector* and returns a single scalar index rather than an `m`-vector of per-column indices — so for `m>1` query points, `idx` silently collapses to length 1, the embedding-selection `mask` logic touches only one of the `m` points, and the rest of `wp`/`fp` are left at their zero-initialized default with no error or warning. Confirmed directly against MATLAB (a single query point through `crmap` matches `crmap0` exactly; batching it with a second point zeros out the second point's result, and 5-point batches error out inside `crembed` entirely). This is a real bug in the original 1998 source for a genuinely degenerate case, not a generator or C++-port issue, and the C++ port's per-column `idx` selection does not (and should not) reproduce it. Sidestepped in `gen_crdiskmap_private.m` by using a single query point for `n3==1` polygons and reserving multi-point batches for `n3>1` polygons, where the bug cannot trigger — consistent with the project's standing preference for fixing/avoiding poorly-conditioned generator inputs over either "fixing" 1998 MATLAB source or papering over it in the port.

### 9.4 `deparam`'s parameter problem can have a near-degenerate Jacobian, making MATLAB's own root choice context-dependent

While porting `deparam` (`@extermap/private/deparam.m`), `gen_extermap_private.m`'s original unit-square test case turned out to be doubly bad: it's both rotationally *and* reflectively symmetric (all four `beta` equal), so its parameter-problem Jacobian is exactly singular at the solution, and *consecutive* calls to `deparam(w,beta,...)` with bit-identical inputs (confirmed via hex dumps of every argument) reproducibly returned different — but equally valid — prevertex orderings depending on whether the call was made directly/via `feval` versus from within `extermap(p,opt)`'s constructor. Replacing the square with an asymmetric quadrilateral did **not** fully fix this: the same direct-vs-constructor split persisted (and was also present, previously unnoticed, in the existing pentagon test case), just manifesting as a different-but-valid root rather than an obviously-symmetric one. Calling `deparam` directly (via `extermap.private_`, matching how `dparam`/`hpparam`'s goldens already call their param functions directly rather than extracting from a constructed map object) is, on its own, fully deterministic and reproducible (10+ repeated calls, identical result every time) — only the constructor route is context-sensitive. Root-caused to MATLAB's own `deparam`/`nesolve` internals (the parameter problem's Jacobian is close enough to degenerate that low-level floating-point execution-path differences between interpreted/`feval` and compiled call contexts select a different valid root), not a generator or C++-port bug. Fixed by changing `gen_extermap_private.m` to compute `deparam`'s literal input the same way `extermap.m`'s constructor does (`flipud(vertex(p))`, `1-flipud(angle(p))`, `scfix('de',...)`) and call `deparam` directly via `extermap.private_` for `z`/`c`/`qdat`, rather than constructing `extermap(p,opt)` and extracting `prevertex`/`constant`/`qdata` from the result; the downstream `dequad`/`dederiv`/`demap`/`deinvmap` cases then un-negate/un-flip (`w = flipud(w_de); beta = 1-flipud(beta_de);`) to recover the normal convention they expect, mirroring `poly = polygon(flipud(w),1-flipud(beta))` in `extermap.m`.

Separately (also caught here): the golden generator had been extracting `deparam`'s "beta" *input* as `angle(m.polygon)-1` (the normal convention used by the downstream functions), but `extermap.m` actually calls `deparam` itself with `beta = 1-flipud(angle(poly))` — the *negation* of that, needed so `sccheck`'s orientation-sum check (`sum(beta)` must equal `+2`, not `-2`, for the clockwise-traversed `'de'` type) passes. The previously-recorded golden `beta` could not actually reproduce the recorded `z`/`c`/`qdat` via a direct `deparam` call; this is fixed by the same generator rewrite above.

### 9.5 `rparam`'s golden generator had the same input/output mismatch as `deparam`'s — fixed the same way

`gen_rectmap_private.m`'s `cases_rparam` paired `w`/`beta` extracted from `m.polygon` (which `@rectmap/rectmap.m`'s constructor stores *after* applying `scfix('r', w, beta, corner)`, i.e. the post-`scfix` values) with the *original*, pre-`scfix` `corners` local variable, never updated to match. This is silent and harmless whenever `scfix` happens not to renumber anything (the `hex6`/`corners 1:4` case, where `corner(1)` is already `1`) — exactly why it went unnoticed — but breaks whenever `scfix` actually renumbers (the `L-shape`/`corners [2,4,5,1]` case, where `scfix` maps the corners to `[1,3,4,6]`). Confirmed directly: calling `rparam` with the recorded `w`/`beta` but the *original* `[2,4,5,1]` corners converges to a different, internally-consistent-but-wrong rectangle (`L≈0.40` vs the recorded `L≈2.12`) — not a near-degenerate-root issue like §9.4, just a plain mismatched pair of inputs. Fixed in `gen_rectmap_private.m` by recomputing `[~,~,corners_post] = scfix('r', vertex(p), angle(p)-1, corners)` and storing `corners_post` (matching the precedent set by §9.4's `deparam` fix and the `dparam`/`hpparam` generators, which call `scfix` themselves rather than relying on a constructed map object's possibly-divergent internal state).

The C++ `rparam` port additionally needed `rptrnsfm.m`'s Trefethen-style "long edge / short edge" transform (long edges filled by direct cumulative-sum-of-exponentials; short edges via a log/exp mountain-pass transform with a near-zero-argument fallback) to convert the unconstrained nesolve variable `y` into strip-image prevertices, and reuses `rectAugQdat`/`rectStripEnds` (built for `rmap`/`rderiv` in Phase 9) for the `stquad`/`stquadh`-based residual integrals. After `nesolve` converges on the strip, a second, separate fixed-point Newton iteration (using `r2strip`'s forward map and derivative) refines the strip-image prevertices onto the actual rectangle boundary, with per-side boundary-clamping logic (real part capped on top/bottom edges, imaginary part capped on left/right edges) to keep iterates from leaving the rectangle. The L-shape golden case needs a slightly relaxed comparison tolerance (5× the case's recorded `tol`) in the C++ test, attributed to ordinary cross-implementation floating-point path divergence compounding across the two sequential iterative solves (`nesolve` then the rectangle-boundary Newton refinement) rather than any structural discrepancy — the `hex6` case (no compounding needed, converges in fewer effective degrees of freedom) matches to the full recorded tolerance.

### 9.6 `crparam` needed a from-scratch triangulation pipeline; one branch is deliberately unported

`crparam` (`@crdiskmap/private/crparam.m`) is the only `XXparam` function that doesn't just solve a nonlinear residual against precomputed quadrature — it first has to *build* the disk-formulation's combinatorial structure (the quadrilateral graph `Q`) from the polygon geometry. This needed five new files ported essentially independently before `crparam` itself could be written:

- **`crtriang`** (`@crdiskmap/private/crtriang.m`): a stack-based recursive ear-clipping triangulator. Picks the sharpest (most acute) vertex, tries the trivial triangle at that vertex, and falls back to a "visibility" diagonal search (handling borderline on-edge points via `crpsdist`, including a slit/crack special case that nudges the test point inward before re-checking `isinpoly`) when other vertices fall inside the trial triangle. Ported using a `std::vector`-based LIFO stack (rather than MATLAB's preallocated dense 0/1-marker-row stack with manual growth) since the *order* of stack push/pop — not the storage mechanism — is what has to match MATLAB bit-for-bit for the resulting triangulation to be combinatorially identical.
- **`crcdt`** (`@crdiskmap/private/crcdt.m`): constrained-Delaunay edge flipping over the initial triangulation, using the standard "sum of the four angles at the diagonal's endpoints `< pi`" flip criterion.
- **`crqgraph`** (`@crdiskmap/private/crqgraph.m`): converts the triangulation into the quadrilateral graph `Q` (`qlvert`/`qledge`/`adjacent`), including the CCW/CW orientation-fixing checks via `scangle`. `QGraph` (`cpp/include/sctoolbox/qgraph.hpp`) gained an `edge` field for this (previously only `qlvert`/`qledge`/`adjacent`, since `crembed`/`crspread`/`crgather` never needed it).
- **`crossrat`** (`@crdiskmap/private/crossrat.m`): trivial — the target crossratios implied by `w` and `Q`.
- **`craffine`** (`@crdiskmap/private/craffine.m`): a second, independent stack-based graph traversal (over `Q.adjacent`) that propagates a single 2-parameter affine correction outward from one quadrilateral to all others, composing affine maps via least-squares fits (`Eigen::ColPivHouseholderQR` standing in for MATLAB's `\` on a 3-point-vs-2-unknown overdetermined system) at each step. Restricted to the case where `w` is fully known (crparam's only call site), so the "deduce missing polygon side" logic in the original `.m` is not ported.
- **`crfixwc`**: straightforward once `crembed`/`crimap0`/`isinpoly` already existed from earlier phases.

**`nesolve` needed one new capability for `crparam`'s residual:** `crparam.m` calls `sctool.nesolvei`, not `sctool.nesolve` — and diffing the two `.m` files shows the *only* behavioral difference is that `nesolvei`'s very first Jacobian is the identity matrix instead of a finite-difference approximation (subsequent restart Jacobians, if any, still use `nefdjac`, same as `nesolve`). Rather than duplicating the ~600-line solver, `sctoolbox::nesolve` gained an `identityInitialJacobian` defaulted-`false` parameter (`cpp/include/sctoolbox/nesolve.hpp`) that selects exactly this one substitution.

**Deliberately not ported: `crsplit`'s dynamic re-triangulation surgery.** `crsplit.m` (called by `crparam` to preprocess the polygon) has two phases: chopping very sharp corners (turning-angle `< -0.75`), then iteratively detecting and splitting "narrow channel" edges using a geodesic-distance check (`crpsgd`, a nested function). Both phases' *detection* logic is fully ported (`cpp/src/crsplit.cpp`). But if phase 2 actually finds a channel that needs splitting, the original `.m` performs intricate live surgery on the in-progress CDT — inserting new vertices/edges/triangles and patching `edge`/`triedge`/`edgetri`'s cross-references in place (`@crdiskmap/private/crsplit.m` lines ~103–143). Confirmed directly against MATLAB that **neither** existing `crdiskmap` test polygon (the unit square or the irregular quadrilateral) ever triggers a split (`crsplit` returns `w` and `orig` unchanged, i.e. a no-op, for both). Given zero test coverage for that branch and the risk of silently-wrong mesh surgery being worse than no implementation, the C++ port throws `std::runtime_error` with a descriptive message if a split is ever detected as necessary, rather than guessing at untested logic. Anyone hitting this in practice would need a polygon with a sufficiently narrow concave channel — neither of the two crdiskmap test fixtures qualifies, and no new fixture requiring it was introduced.

With this, the full `cases_crparam` golden (`crdiskmap_private.gold`) — both the unit square (trivial `cr=1`) and the irregular quadrilateral (`cr=[1.5026, 3.4842]`) — matches exactly on the first attempt once each component above was individually cross-checked against MATLAB (`crtriang`/`crcdt`/`crqgraph`'s intermediate `Q.qlvert`/`Q.qledge`/`Q.adjacent` were diffed directly against MATLAB for both fixtures before assembling `crparam` itself, catching the row-order-within-`triedge`-columns non-issue early: `crtriang`'s `triedge` can come out as a harmless cyclic rotation of MATLAB's per-triangle edge ordering, since `crqgraph` always re-rotates each triangle's edge list to put the diagonal first before using it, making the absolute starting offset unobservable downstream).

## 10. Phase 11 — Map Classes

`DiskMap`/`HplMap`/`ExterMap`/`StripMap`/`RectMap`/`CrDiskMap` (`cpp/include/sctoolbox/*map.hpp`) wrap the already-ported `XXparam`/`XXmap`/`XXinvmap`/`XXderiv` functions with the four user-facing methods every MATLAB map class exposes: `eval` (forward map, with each class's own validity filter — `abs(zp)<=1+eps` for disk/exterior/cross-ratio-disk, `imag(zp)>-eps` for half-plane, `0<=imag(zp)<=1` for strip, `isinpoly` against the corner rectangle for rectmap, `isinpoly` against the polygon for cross-ratio-disk's `evalinv`), `evalinv` (inverse map via ODE continuation + Newton polish, using `accuracy()` as the default tolerance exactly as the `.m` classes do), `evaldiff` (derivative), and `accuracy` (a self-consistency estimate, memoized after first computation, with an algorithm specific to each domain — disk/half-plane/exterior compare `c*∫f'` between consecutive prevertices to the corresponding vertex difference; strip/rectangle additionally splice in the strip's two infinite ends and a cross-strip check pair; cross-ratio disk instead compares actual vs. target *crossratios* directly, with no integration-vs-difference comparison at all). Each class's constructor restricts to the "polygon in, parameter problem solved with an automatic initial guess" form, matching every `XXparam` port's own restriction — the MATLAB classes' continuation-map, given-prevertex, and tol-override constructor branches are not ported.

`ExterMap` and `RectMap` needed real care beyond mechanical wrapping:

- **`ExterMap`** re-derives the same negated/flipped `(w,beta)` convention as `deparam.m` internally (`w=flipud(vertex(poly))`, `beta=flipud(1-angle(poly))`) in *every* method, not just the constructor — confirmed by direct translation of `@extermap/eval.m`/`evalinv.m`/`evaldiff.m`/`accuracy.m`, all of which literally re-derive this convention from `polygon(M)` on each call rather than caching it. Got this wrong once while writing this class's own golden generator (storing `1-flipud(beta_de)`, i.e. the *alpha* convention, under the field name `beta`) — caught immediately by the resulting "Invalid polygon" exceptions when reconstructing a `Polygon` from it.
- **`RectMap::accuracy()`** is the most involved of the six (mirroring `rparam`'s own complexity): it renumbers via `rcorners`, maps prevertices to the strip via `r2strip`, splits prevertices into bottom/top groups using the renumbered third corner index, adds one cross-strip check pair, and splices in the strip's `+Inf`/`-Inf` ends (reusing `rectStripEnds`/`rectAugQdat` from Phase 9) before integrating with `stquad`.
- Both `gen_rectmap_class.m`'s and (already-existing) `gen_rectmap_private.m`'s golden generators needed the same fix documented in §9.5: pairing post-`scfix` `w`/`beta` (extracted from a constructed map's `m.polygon`) with the *pre*-`scfix` `corners` argument silently works only when `scfix` happens not to renumber (`hex6`), and is wrong whenever it does (`L-shape`) — fixed by recomputing and storing the post-`scfix` corners alongside `w`/`beta` in every case that needs to reconstruct a fresh `RectMap`.

**A real bug was found and fixed in `ode45.cpp` along the way**, while debugging an unrelated `StripMap::evalinv` golden mismatch for one specific query point: when the ODE right-hand side evaluates to NaN during a trial step (e.g. because the adaptive stepper's trial point landed far enough from the true solution to approach a critical point of the conformal map, where the derivative — and hence `1/derivative` in the inverse-map ODE — blows up), the resulting `errNorm` is NaN. The old step-size-adaptation code computed `factor = (errNorm > 0.0) ? shrink-formula : 5.0` — and `errNorm > 0.0` is `false` for NaN, identical to the *legitimate* `errNorm==0` ("step was perfect, take a bigger one") case, so a NaN step caused the integrator to **grow** the step size instead of shrinking it, and (separately) a NaN step could even be silently *accepted* once `h` underflowed below the `1e-14` escape threshold. Fixed by explicitly detecting `std::isnan(errNorm)` and forcing an aggressive shrink (and never accepting that step), plus adding a `kMaxSteps` safety counter so a persistently-NaN RHS can't hang the integrator. This fix did not, on its own, rescue the specific failing `StripMap::evalinv` query point — that turned out to be a genuine ODE-path sensitivity (MATLAB's `ode23` happens to take a different adaptive path that avoids passing near the same critical point; this port's Dormand-Prince `ode45` does not, for that specific starting point) rather than a bug, and was resolved the same way as the analogous `deinvmap` issue in §9.1: by choosing better-conditioned query points in the golden generator. The `ode45.cpp` fix is still a correctness improvement worth keeping regardless, since the backwards NaN-handling could in principle have caused silent corruption (an accepted-but-wrong step) rather than a clean NaN propagation in some other, as-yet-unexercised scenario.

## 11. Phase 12 — `annulusmap` (doubly-connected D-SC map)

`AnnulusMap` (`cpp/include/sctoolbox/annulusmap.hpp` + `annulusmap_internal.hpp`) ports the DSCPACK algorithm (Chenglie Hu's 1995 "User's Guide to DSCPACK") used by `@annulusmap` for conformal maps between the canonical annulus `{u<=|w|<=1}` and a doubly connected polygonal region. Restricted, like every other map class in this codebase, to the bounded-outer-polygon constructor form — the MATLAB `'truncate'` option for unbounded outer polygons (`ishape==1` throughout `qinit`/`dscsolv`/`dscfun`) is not ported and throws if reached.

The dependency chain mirrors the `XXparam`/`XXmap`/`XXinvmap` pattern but is built around Jacobi theta-functions instead of Gauss-Jacobi-quadrature-of-a-product-of-powers: `qinit` (Gauss-Jacobi nodes/weights, flat 1D work array `qwork`, kept 0-indexed but with all the loop arithmetic translated 1:1 from MATLAB's 1-indexed formulas rather than re-derived, to avoid off-by-one mistakes) → `thdata`/`wtheta` (theta-function series/closed-form evaluation, branching on the inner radius `u`: series sum for `u<0.63`, a `cosh`-based closed form for `0.63<=u<0.94`, and a further-simplified closed form for `u>=0.94`) → `wprod` (the D-SC integrand, a product of theta-function ratios) → `wqsum`/`wquad1`/`wquad` (compound one-sided Gauss-Jacobi quadrature along line segments or circular arcs, with adaptive subdivision driven by distance-to-nearest-prevertex) → `xwtran` (unconstrained-parameter-vector → actual D-SC parameters `u,c,w0,w1,phi0,phi1`) → `dscfun` (the nesolve residual) → `dscsolv` (the parameter-problem driver) → `zdsc`/`wdsc` (forward/inverse map, the latter via Euler-scheme initial guess + Newton polish + recursive retry with a perturbed reference vertex on failure, exactly matching `wdsc.m`).

**Two real bugs were found and fixed while debugging this phase, both via a "does this match MATLAB's literal closure capture / argument list" line-by-line diff rather than re-deriving the math from scratch:**

1. **`dscfun`'s fixed prevertex was silently zero instead of 1.** MATLAB's `dscfun.m` is a nested closure that captures `w0`/`phi0` from `dscsolv.m`'s enclosing scope, where `w0(M)=1`/`phi0(M)=0` was set *once* before the nesolve loop began and never touched again by `xwtran` (which only ever writes indices `1..M-1`). The direct C++ translation allocated a fresh zero-initialized `w0`/`phi0` *inside* `dscfun` on every call, leaving `w0(M-1)` at `0` instead of `1`. Since `w0(M-1)` is used as a quadrature endpoint and (via `u*w0` inside `wprod`) as a divisor, this didn't crash — it produced a finite but wrong residual, which sent `nesolve` first into a 200+ call non-convergent spiral and (before a debug call-counter was added to find that) looked indistinguishable from an honest infinite loop. Fixed by setting `w0(M-1)=1`, `phi0(M-1)=0` explicitly at the top of the C++ `dscfun`, matching the closure-captured constant.
2. **Wrong quadrature radius for the inner-polygon equations.** `dscfun.m`'s `win1` (the inner-polygon rotation-fixing equation) and its `N-1` inner-polygon side-length conditions both call `wquad(...)` with the `radius` argument (9th positional parameter) set to `u` — the inner prevertex circle's radius — since the integration path runs around the inner circle. The initial C++ translation passed `1.0` (the outer circle's radius) at both call sites, a transcription slip from the visually-similar `win3`/outer-side-length calls a few lines below (which correctly use radius `1.0`). This was diagnosed by writing a temporary `debugDscfun` export and a parallel MATLAB script that reconstructed `dscsolv`'s exact initial-guess `x` and called `annulusmap.private_('dscfun', x, fdat)` directly: the two 10-element residual vectors matched exactly in every entry *except* the `win1`/inner-side-length block (entries 1–5), pinpointing the bug immediately. After both fixes, `nesolve` converges in the same shape as MATLAB (symmetric `u=0.5412...`, `w0=[i,-1,-i,1]`, `phi0` spaced exactly `pi/2` apart for the symmetric square-annulus test fixture) and the forward/inverse map goldens match to the recorded tolerances.

`gen_annulus_private.m`'s existing skeleton (from an earlier session) had its own latent bug, caught while writing the C++ test: it computed the eval/evalinv query-point radii as `exp(-u) + frac*(1-exp(-u))`, treating `u` as if it were a *log*-radius — but `u` is the inner radius of the canonical annulus directly (`|w1(k)|==u` by construction in `xwtran.m`). Fixed to `u + frac*(1-u)`. The generator was also extended to store the outer/inner polygon's raw vertex/angle arrays (not just the solved `u,c,w0,...` fields, which `annulusmap`'s restrictive `subsref` blocks from being read outside `golden_fields.m` anyway) so the C++ test can reconstruct an `AnnulusMap` end-to-end via its public constructor — exercising the *whole* pipeline (`qinit`→`dscsolv`→`eval`/`evalinv`) per case, the same way `gen_crdiskmap_class.m` etc. do for the simply-connected map classes, rather than only the lower-level private functions.

`wquad1.m`'s line-segment (`linearc==0`) continuation branch contains a self-indexing typo (`d(d(d ~= 0))`, using distance *values* as array *subscripts*) in its `while` loop body that would throw a MATLAB indexing error if ever executed; like `crsplit`'s deferred mesh-surgery branch (§9.6), this is not ported — the C++ `wquad1` throws `std::runtime_error` if that continuation is ever needed, since no exercised test case hits it and guessing at the intended fix was judged worse than failing loudly. The circular-arc (`linearc==1`) continuation loop, which *is* exercised, has a `kMaxSteps`-style safety counter (10000 iterations) added defensively, matching the precedent set by the `ode45.cpp` fix in Phase 11.

## 12. Phase 13 — `moebius` and `composite`

`Moebius` (`cpp/include/sctoolbox/moebius.hpp`) ports `@moebius` in full: the four-coefficient and three-point constructors, `eval`, `diff`, `inv`, `normal`, the `M1(M2)` composition form of `subsref.m`, and the scalar arithmetic operators (`plus`, `minus`, `mtimes`, `mrdivide` in both operand orders, `uminus`, `uplus`). Only the pretty-printers (`char`/`disp`/`display`) are skipped. The pre-existing `moebius3` (`cpp/src/moebius3.cpp`) is left alone: it is the all-finite fast path of the same three-point formula, called on `crspread`/`crgather` hot paths, and folding it into the new class would have churned already-passing Phase 10 goldens for no gain.

The three-point constructor's infinity handling is the only intricate part, and it is ported branch-for-branch rather than re-derived: the `rem((j-1:j+1)+2,3)+1` renumbering that rotates the infinite entry into the middle slot, the separate `w`-infinite and `z`-infinite paths, the nested "move `Inf` to the beginning of `z`" swap, and the `isnan(A(1))` sentinel that MATLAB uses to detect whether the special-case branch already filled in the coefficients (a `bool haveA` in C++). The renumbered `z`/`w` are stored as `source`/`image`, since `inv()` swaps them. `gen_moebius.m` covers all thirteen reachable combinations of which slot (if any) of `z` and `w` holds an infinity.

**`eval.m` and `diff.m` disagree about infinity, and the goldens record that.** `eval.m` deliberately funnels every degenerate result to `Inf`: infinite inputs map to `c2/c4`, a denominator with `|den| < 3*eps` is forced to `NaN`, and a trailing `f(isnan(f)) = Inf` sweep converts every `NaN` — however it arose — into `Inf`. `diff.m` (added in 2007, nine years after the rest of the class) has no infinity handling at all, so at `z = Inf` it evaluates `0*Inf` and returns a genuine `NaN`. This is not a porting artifact; it is what MATLAB returns, and `moebius.gold` stores `fp = Inf` alongside `dp = NaN` for the same input point. `test_moebius.cpp` therefore compares non-finite entries in two modes: `kInfExact` for `eval` outputs (MATLAB never returns `NaN` there, so `Inf` must match `Inf`), and `kAnyNonFinite` for `diff` outputs. The looser mode for `diff` is deliberate — C++ agrees with MATLAB and yields `NaN` on this platform, but the standard explicitly permits a complex multiplication to recover an infinity where IEEE arithmetic would produce `NaN`, so asserting which flavor of non-finite comes out would be testing the standard library rather than this port.

**`Composite` needs explicit type erasure, and that is what makes `inv()` work.** MATLAB's `composite.m` is duck-typed: it accepts `moebius`, any of the six SC map classes, `scmapinv`, or an `inline`, stores them in a cell array, and dispatches through `feval`. C++ has no equivalent, so `Composite::Member` holds a forward callable *plus* its inverse callable (empty when there is none). That pairing is not an implementation convenience — it is the direct analogue of what MATLAB does: `inv(composite)` reverses the member list and calls `inv` on each member, and `inv` of an SC map yields an `scmapinv` whose `eval` is the original map's `evalinv`. Reversing the list and swapping `forward`/`inverse` per member reproduces this exactly. Members are built with the factories `Composite::member` (templated over any class with `eval`/`evalinv`, plus a `Moebius` overload), `Composite::inverseMember` (the `scmapinv` case), and `Composite::function` (the `inline` case — usable but not invertible, so `inverse()` throws `std::runtime_error`, mirroring MATLAB's "Can't invert INLINE maps."). Maps are copied into the member closures, so a composite outlives the objects it was built from. `append(const Composite&)` splices in another composite's members, matching the flattening MATLAB's constructor performs.

`gen_composite.m` pins three member orderings per polygon (`map_then_mob`, `mob_then_map`, `mob_map_mob`), recording the ordering as a `STR` golden field so the C++ test can rebuild the same chain. The Moebius map placed *before* a `diskmap` is a Blaschke factor `(z-a)/(1-conj(a)z)`, i.e. a disk automorphism, since anything else would feed the disk map points outside its domain.

**`generateGoldens.m`'s `ALL` list was missing the six `_class` groups.** `gen_diskmap_class.m` and its five siblings existed and their `.mat`/`.gold` outputs were committed, but the group names were never added to `ALL`, so `generateGoldens()` with no arguments silently skipped them and `generateGoldens('diskmap_class')` errored as an unknown group. Fixed while adding `moebius`/`composite`; `ALL` is now the complete list of 21 groups, and a full `generateGoldens()` followed by `exportGoldensText(ALL)` reproduces every committed `.gold` file byte-for-byte.
