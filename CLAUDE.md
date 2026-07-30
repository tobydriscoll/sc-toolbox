# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

**Run all tests** (from MATLAB, with the repo root on the path):
```matlab
result = sctool.runTests()
```

**Run a single map type's tests:**
```matlab
result = sctool.runTests('disk')       % testDisk
result = sctool.runTests('ext')        % testExterior
result = sctool.runTests('hpl')        % testHalfplane
result = sctool.runTests('strip')      % testStrip
result = sctool.runTests('rect')       % testRectangle
result = sctool.runTests('crdisk')     % testCRDisk
result = sctool.runTests('ann')        % testAnnulus
```

**Regenerate fixture golden values** (only needed when polygons or solver logic changes):
```matlab
generateFixtures()               % all types
generateFixtures('disk')         % one type
generateFixtures({'disk','ann'}) % subset
```
Fixtures are stored in `tests/fixtures.mat`. Tests skip (not fail) if this file is absent.

**Add the toolbox to the MATLAB path** (must be done once per session):
```matlab
addpath(genpath('/path/to/sc-toolbox'))
```

### C++ port

Prerequisites: a C++17 compiler and `cmake`. Eigen and Catch2 are fetched by CMake (pinned by version + hash), so no other installs are needed. MATLAB is **not** needed to build or test the C++ code — the goldens are committed as text.

**Build and run the C++ tests:**
```sh
cmake -S cpp -B cpp/build
cmake --build cpp/build
./cpp/build/sctoolbox_tests     # or: ctest --test-dir cpp/build
```

**Regenerate the C++ golden values** (MATLAB; only when a `.m` source or a generator changes):
```matlab
generateGoldens()      % -> tests/cpp/goldens/*.mat
exportGoldensText()    % -> tests/cpp/goldens_text/*.gold  (what the C++ tests read)
```

Status, phase-by-phase progress, and the remaining work live in `CPP_PLAN.md` §0.

## Architecture

### Class hierarchy

`scmap` is the abstract base class (polygon + options). Seven map types inherit from it:

| Class | Domain | Target |
|---|---|---|
| `diskmap` | unit disk | polygon interior |
| `hplmap` | upper half-plane | polygon interior |
| `extermap` | disk exterior | polygon exterior |
| `stripmap` | infinite strip | polygon interior |
| `rectmap` | rectangle | generalized quadrilateral |
| `crdiskmap` | disk (cross-ratio formulation) | polygon interior |
| `crrectmap` | rectilinear polygon | polygon (cross-ratio) |

`annulusmap` is a standalone class (not inheriting `scmap`) for doubly-connected regions (outer + inner polygon boundary). `moebius` and `composite` are also standalone utility classes.

### Per-class private function pattern

Each `@maptype/` folder has a `private/` subfolder with the low-level numerical kernel:

- `XXparam.m` — solves the parameter problem (calls `sctool.nesolve`)
- `XXmap.m` — forward SC integral evaluation
- `XXderiv.m` — derivative (product formula `c * ∏(z - z_k)^β_k`)
- `XXinvmap.m` — inverse map (ODE continuation + Newton polish)
- `XXimapfun.m` — ODE right-hand side for `XXinvmap`
- `XXquad.m` — Gauss-Jacobi quadrature on the preimage domain

`rectmap` additionally contains `ellipjc.m` (complex Jacobi elliptic functions via Landen transformation) and `ellipkkp.m` (K and K′ via AGM), since the rectangle map requires elliptic integrals.

### Quadrature system

All SC integrals use **Gauss-Jacobi quadrature** tailored to the algebraic singularity at each prevertex.

- `sctool.gaussj(n, alf, bet)` — nodes and weights via Lanczos iteration + symmetric tridiagonal `eig`
- `sctool.scqdata(beta, nqpts)` — builds the `qdat` matrix (`[nodes, weights]`, one column pair per vertex) consumed by every `XXquad.m`
- Each `XXquad.m` performs adaptive subdivision: if a singularity lies within half the interval width of the left endpoint, the interval is split and the first piece uses Gauss-Jacobi; remaining pieces use regular Gauss quadrature

### Nonlinear solver

`sctool.nesolve` is a self-contained Newton solver (Dennis & Schnabel 1983, Algorithm D6.1.3). It supports two globalization strategies selected by `details(2)`:
- `1` — line search with cubic backtracking (`nelnsrch`)
- `2` — hook step / trust-region (`nehook`, `netrust`)

Support files: `nefdjac` (finite-difference Jacobian), `neqrdcmp` (Householder QR), `nechdcmp` (perturbed Cholesky), `nemodel` (affine model), `nebroyuf` (Broyden rank-1 update), `nestop` (termination).

The annulus parameter problem (`dscsolv`) also calls `nesolve`.

### Polygon normalization

`scfix(type, w, beta)` is called before every parameter problem. It enforces domain-specific ordering constraints (e.g., for disk/half-plane: vertices 1, 2, n−1 finite; β(n) ∉ {0,1}) by cyclically renumbering vertices and, if necessary, inserting zero-turn vertices via `scaddvtx`.

### Inverse maps

All inverse maps use **ODE continuation followed by Newton refinement**:
1. Integrate the inverse ODE (using `ode23` for hpl/strip/exterior/rect/crdisk, `ode113` for disk) from a known base point to the target
2. Polish each result with Newton steps using `XXmap` as the residual and `XXderiv` as the Jacobian

### Options

`sctool.scmapopt` (or `scmapopt` on the `@scmap` class) constructs the options struct. Key fields: `Tolerance` (parameter problem residual), `Trace` (verbosity), `InitialGuess`, `Method` (nesolve vs fsolve).

## Tests

Tests live in `tests/` as `matlab.unittest.TestCase` subclasses. Each class loads `tests/fixtures.mat` in `TestClassSetup` and skips all tests if the file is absent. Test methods verify forward eval, inverse eval, derivative, accuracy estimate, and finite-difference derivative consistency. Tolerances are `fixture.tol * 50` (default `tol = 1e-12`, so effective tolerance ≈ 5e-11).

The C++ translation plan and function-level golden-value test infrastructure are documented in `CPP_PLAN.md`.
