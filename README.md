# EndoBeams

`EndoBeams.jl` is a Julia finite element package for beam-to-surface contact problems. The package is based on a 3D FE corotational formulation for frictional contact dynamics of beams where the target surface, supposed rigid, is described implicitly using a Signed Distance Field (SDF), predefined in a volumetric grid.

----------------------------

## Basic usage
First, add the Endobeams.jl package with `using Pkg, Pkg.add(url="https://github.com/beatricebisighiniEMSE/EndoBeams.jl")`.
Run one of the examples present in the "examples" folder:
- `angle.jl`: Right-angle cantilever beam subject to out-of-plane loading;
- `ring.jl`: Impact of a ring against a rigid surface;
- `net.jl`: Dropping a net on a rigid sphere;
- `stent.jl`: Deployment of braided stent.

----------------------------
## Branches
- `main`: The latest package version;
- `original-version-article`: Old version - the package version from [1];
- `FEM-stent-deployment`: Old version - stick-slip friction + all the utils to perform patient-specific simualation of braided stent deployment (one case provided);
- `beam-to-beam`: Old version - first attempts with beam-to-beam contact formulation.

----------------------------

## References
[1] Aguirre M, Avril S. 2020. An implicit 3D corotational formulation for frictional contact dynamics of beams against rigid surfaces using discrete signed distance fields. Comput Methods Appl Mech Eng. 371:113275.
[2] Bisighini, B., Aguirre, M., Pierrat, B., & Perrin, D. (2022). Advances in Engineering Software EndoBeams . jl : A Julia finite element package for beam-to-surface contact problems in cardiovascular mechanics. 171(July). https://doi.org/10.1016/j.advengsoft.2022.103173
