module EndoBeams

include("packages.jl")
include("constants.jl")
include("structures.jl")
include("type_node.jl")
include("type_beam.jl")
include("solver.jl")
include("contact.jl")
include("constraints.jl")
include("utils_solver.jl")
include("sparse_matrices.jl")
include("fem.jl")
include("params.jl")
include("IO.jl")

export constructor_nodes, constructor_beams, constructor_constraints, Params, ParamsTest, solver!, Material, Geometry, SimulationParameters, ExternalForces, BoundaryConditions, Configuration, Vec3, Vec2, Mat33, Discrete_SDF, Plane_z_SDF, Plane_y_SDF, Sphere_SDF, Cylinder_SDF


include("precompile/precompiles.jl")
_precompile_()

end
