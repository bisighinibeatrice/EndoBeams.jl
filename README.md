# EndoBeams

`EndoBeams.jl` is a Julia finite element package for beam-to-surface contact problems. The package is based on a 3D FE corotational formulation for frictional contact dynamics of beams where the target surface, supposed rigid, is described implicitly using a Signed Distance Field (SDF), predefined in a volumetric grid.

## 🚀 Features

- Corotational beam elements for modeling quasi-inextensible wire structures with large displacements
- Implicit surface description via signed distance fields (SDF) on a regular voxel grid
- Frictional contact handling with penalty and regularized friction models
- Highly optimized Julia implementation (faster than MATLAB, comparable or better than Abaqus for specific applications)
- Built-in examples for cantilever bending, impact, net drop, and braided stent deployment

## 📦 Installation

Install using Julia's package manager:

```
using Pkg
Pkg.add(url="https://github.com/bisighinibeatrice/EndoBeams.jl")
```

## ▶️ Examples
Example simulations are located in the `examples/` directory:

- `angle.jl` — Cantilever beam under vertical loading
- `ring.jl` — Ring dropping onto a rigid surface
- `net.jl` — Net dropped onto a rigid sphere
- `stent.jl` — Braided stent deployed into a rigid cylindrical vessel

Each script includes:

- Mesh reading or writing 
- Node and element construction
- Contact and solver setup
- Simulation loop
- Result export to `.vtk` files (for ParaView)

## 📁 Repository Structure

```
EndoBeams.jl/
├── src/             # Core finite element and contact modules
├── examples/        # Ready-to-run simulations
├── test/            # Unit tests
├── Project.toml     # Julia environment declaration
├── Manifest.toml    # Package dependency snapshot
└── LICENSE          # MIT License
```

## 🔀 Branches

- `master` — current maintained version
- `original-version-article` — original version matching 2022 publication
- `stent_deployment` — includes codes to performe stent-deployment simulation (branch of `master`)
- `beam-to-beam` — includes codes to model beam-to-beam contact (branch of `original-version-article`)

## 📚 References
[1] Aguirre M, Avril S. 2020. An implicit 3D corotational formulation for frictional contact dynamics of beams against rigid surfaces using discrete signed distance fields. Comput Methods Appl Mech Eng. 371:113275. https://doi.org/10.1016/j.cma.2020.113275

[2] Bisighini, B., Aguirre, M., Pierrat, B., & Perrin, D. (2022). Advances in Engineering Software EndoBeams . jl : A Julia finite element package for beam-to-surface contact problems in cardiovascular mechanics. 171(July). https://doi.org/10.1016/j.advengsoft.2022.103173

## 🤝 Contact
beatrice.bisighini@emse.fr
beatricebisighini@gmail.com