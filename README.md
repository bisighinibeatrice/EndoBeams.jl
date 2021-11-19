# EndoBeams

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](http://pierrat.gitlab.emse.fr/EndoBeams.jl/dev)
[![Build Status](https://gitlab.emse.fr/pierrat/EndoBeams.jl/badges/master/pipeline.svg)](https://gitlab.emse.fr/pierrat/EndoBeams.jl/pipelines)
[![Coverage](https://gitlab.emse.fr/pierrat/EndoBeams.jl/badges/master/coverage.svg)](https://gitlab.emse.fr/pierrat/EndoBeams.jl/commits/master)

 `EndoBeams.jl` is a Julia finite element package for beam-to-surface contact problems. The package is based on a 3D FE corotational formulation for frictional contact dynamics of beams where the target surface, supposed rigid, is described implicitly using a Signed Distance Field (SDF), predefined in a volumetric grid.

-----------------------------

## Basic usage
First, add the Endobeams.jl package with `using Pkg, Pkg.add(url="https://gitlab.emse.fr/pierrat/EndoBeams.jl")`.
Run one of the examples present in the "examples" folder:
- `angle.jl`: Right-angle cantilever beam subject to out-of-plane loading;
- `ring.jl`: Impact of a ring against a rigid surface;
- `net.jl`: Dropping a net on a rigid sphere;
- `stent.jl`: Deployment of braided stent.

## References
Aguirre M, Avril S. 2020. An implicit 3D corotational formulation for frictional contact dynamics of beams against rigid surfaces using discrete signed distance fields. Comput Methods Appl Mech Eng. 371:113275.