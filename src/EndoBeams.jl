module EndoBeams

include("packages.jl")
include("constants.jl")
include("type_node.jl")
include("type_beam.jl")
include("constraints.jl")
include("contact.jl")
include("structures.jl")
include("linearsolvers.jl")
include("solver.jl")
include("utils_solver.jl")
include("sparse_matrices.jl")
include("fem.jl")
include("params.jl")
include("IO.jl")

export build_nodes, build_beams, build_constraints, Params, solver!, Material, Geometry, ContactParameters, BeamProperties, ExternalForces, BoundaryConditions, Configuration, Vec3, Vec2, Mat33, Discrete_SDF, Plane_z_SDF, Plane_y_SDF, Sphere_SDF, Cylinder_SDF, ID3, get_Rₑ⁰


include("precompile/precompiles.jl")
_precompile_()

end
