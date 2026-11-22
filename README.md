# hdg-maxwell-cylinder

A minimal, educational C++ demo of a hybridizable discontinuous Galerkin (HDG) method for time-harmonic Maxwell scattering by a 2D circular cylinder (TM mode).

> Demo / skeleton only. Simplified P1 / P0 HDG, small polar mesh, dense solve. No results are claimed as validated. Do not use as a production solver or as a reference solution.

## Problem

2D TM mode, time dependence `exp(-i ω t)`:

```
-∇² Ez - k² Ez = 0   in the annulus a < r < R
Ez = 0               on the inner cylinder r = a  (PEC)
∂Ez/∂n + i k Ez = 2 i k Ez_inc   on r = R  (first-order absorbing / Robin)
Ez_inc = exp(i k x),  k = 2π f / c
```

Default parameters in `src/main.cpp`:

- `a = 1.0`, `R = 3.0`
- `f = 1.0 GHz`, `c = 3e8 m/s`
- radial layers `Nr = 4`, angular divisions `Nt = 16`
- stabilization `τ = k` (real, for the demo; try `τ = i k` in code)

## Method (simplified)

- Polar structured mesh, each quad split into 2 triangles.
- Per element: P1 `Ez`, P0 trace `λ` on each edge.
- Local matrix `A = K - k² M + τ E`, coupling `B` from `τ` on edges, `C = τ |edge|`.
- Static condensation (Schur complement) eliminates the interior `Ez`:
  `S = C - Bᵀ A⁻¹ B`, solve the global trace system for `λ`, then recover `Ez` per element.
- Dense complex Gaussian elimination. Small mesh only — this is intentional, to keep the demo readable.

This is the single-trace baseline. A double-trace HDG (both tangential `E` and `H` traces as hybrid unknowns) is a natural extension, noted as a TODO in the code.

## Build

```bash
cmake -S . -B build
cmake --build build -j
./build/hdg_cylinder
```

Or without CMake:

```bash
g++ -std=c++17 -O2 -Iinclude src/main.cpp -o hdg_cylinder
./hdg_cylinder
```

Output:

- console: mesh size, trace unknowns, residual, sample `Ez`
- `field.csv`: `x,y,Ez_real,Ez_imag,|Ez|` at element centroids

## Layout

```
CMakeLists.txt
include/cylinder_mesh.h   polar mesh for the annulus
include/hdg_tm.h          local HDG matrices, condensation, dense solver
src/main.cpp              parameters, assembly, output
```

## License

MIT. See `LICENSE`.
