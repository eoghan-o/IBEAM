# IBEAM

This is the repository for the Incidence-Based Edge-Averaged Method (IBEAM). It is MATLAB code for solving convection–diffusion problems on general polygonal (2D) and polyhedral (3D) meshes. The method uses the mesh's edge–vertex incidence structure to build each element's stiffness matrix, with Bernoulli-function (exponential) edge averaging of the convection term. Wachspress basis functions supply the weights, and the pieces combine as `A_T = Gᵀ W D_α G_b`.

The repository includes 2D and 3D solvers, mesh generators, error and convergence tools, and ParaView export for the 3D results.

## The Problem

The drivers solve the conservative-form convection–diffusion equation

```
-∇ · (α∇u + βu) = f
```

on the unit square (2D) or unit cube (3D), with Dirichlet boundary data taken from a manufactured solution:

- Diffusion coefficient `α` (tested at 1, 0.1, 0.01, and 0.001)
- Constant convection vector `β`, called `b` in the code (`[0; 1]` in 2D, `[0; 0; 1]` in 3D)
- Exact solution `u = sin(πx) sin(πy)` (2D) or `sin(πx) sin(πy) sin(πz)` (3D)
- Source term `f` computed from the exact solution

Errors are measured in the L² norm and the H¹ semi-norm. Smaller `α` makes the problem more convection dominated.

## How It Works

For each element `T`, `ibeam.m` (2D) and `ibeam3d.m` (3D) assemble a local matrix as follows:

1. Build the element's edge list, including **virtual edges** between non-adjacent vertices.
2. Build `G`, the edge–vertex difference matrix.
3. Build `G_b`, which averages across each edge using the Bernoulli function of the local Péclet number `p_e = (β · e) / α`, where `e` is the edge vector.
4. Build the diagonal weight matrix `W`, using Wachspress basis gradients integrated by quadrature over the element.
5. Form `A_T = α · Gᵀ W G_b` and add it into the global sparse matrix.

Boundary values are then imposed and the interior system is solved directly.

## Requirements

- MATLAB. The code uses `sparse`, `convhull`, and `evalc`, with no special toolboxes noted.
- **[PolyMesher](https://paulino.princeton.edu/software)** for the 2D polygonal meshes. It is not included in this repo. Download it and place it in `2d-IBEAM/PolyMesher/` (that path is gitignored).
- [ParaView](https://www.paraview.org/) to view the 3D `.vtu` output (optional).

## Quick Start

### 2D

```matlab
cd 2d-IBEAM
ibeam_2d
```

Edit the parameters at the top of `ibeam_2d.m`:

| Setting | Meaning |
|---|---|
| `k_values` | Mesh refinement levels (`3:9`) |
| `a_values` | Diffusion coefficients (`(1/10).^(0:3)`) |
| `b` | Convection vector |
| `nquad` | Quadrature order |
| `use_poly` | `true` for PolyMesher polygons, `false` for a structured quad mesh |
| `perturb` | Mesh perturbation amount for the structured mesh (`0` gives a regular mesh) |

The script plots the solution, the mesh, and L² and H¹ convergence for each `α`, then saves the results to `polygons.mat`.

### 3D

```matlab
cd 3d-IBEAM
main
```

This runs a refinement study on hexahedral meshes (`generate_mesh3d_hex`) with `k_values = 2:6` and writes `.vtu` files to `3d-IBEAM/Paraview/` named like `ibeam3d_k4_a0.1.vtu`. A tetrahedral generator (`generate_mesh3d`) is also included.

### Convergence Tables

`table_maker.m` reads saved `.mat` result files and generates LaTeX tables, which are written to `convergence_tables.tex`.

## Repository Layout

```
IBEAM/
├── 2d-IBEAM/
│   ├── ibeam_2d.m            # 2D driver: refinement study and plots
│   ├── ibeam.m               # 2D assembly and solve
│   ├── wachspress2d.m        # Wachspress basis functions and gradients
│   ├── assemble_weight_matrix.m
│   ├── polyintegrate.m       # Quadrature over polygons
│   ├── triquad.m             # Triangle quadrature rules
│   ├── generate_mesh.m       # PolyMesher polygonal mesh
│   ├── generate_perturbed_mesh.m
│   ├── generate_hex_mesh.m
│   ├── L2_error.m, H1_semi_error.m
│   ├── polygons.mat          # Saved polygon-mesh results
│   └── Plots/                # Meshes, solutions, convergence plots
├── 3d-IBEAM/
│   ├── main.m                # 3D driver: refinement study and ParaView export
│   ├── ibeam3d.m             # 3D assembly and solve
│   ├── wachspress3d.m        # 3D Wachspress basis functions
│   ├── build_*.m             # Edge, face, graph/normal, G, G_b, and RHS builders
│   ├── bernoulli.m
│   ├── tetquad.m, polyintegrate3d.m
│   ├── generate_mesh3d.m, generate_mesh3d_hex.m
│   ├── L2_error3d.m, H1_semi_error3d.m
│   ├── export_to_paraview.m
│   ├── plot_*.m              # Convergence and wireframe plots
│   ├── testing.m, diagnostics.m
│   └── Paraview/             # Generated .vtu output
├── table_maker.m
└── convergence_tables.tex
```

## Results

Convergence plots for the 2D meshes (squares, perturbed quadrilaterals, and generalized convex polygons) are in `2d-IBEAM/Plots/`. The L² and H¹ errors there are compared against reference slopes of `O(h²)` and `O(h)`.

For general convex polygons at `α = 1`, `convergence_tables.tex` shows L² orders of about 1.7 to 1.9 and H¹ orders of about 0.84 to 0.95.

## Notes and Known Issues

- **ParaView export and hex meshes.** `export_to_paraview.m` writes every cell as a tetrahedron (`VTK_TETRA`), but `main.m` currently uses the 8-node hex mesh. Hex output will not display correctly until the exporter handles hexahedra (VTK type 12) or `main.m` uses the tet mesh.
- **PolyMesher is required** for `use_poly = true` and is not bundled.
- The 3D tables in `convergence_tables.tex` are not filled in yet.

## Acknowledgments

`wachspress2d.m` is credited in its header to M. Floater et al. PolyMesher is by the Paulino research group.
