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
include("utils_fem.jl")
include("params.jl")
include("rotations.jl")
include("IO.jl")

export Params, ParamsTest, constructor_constraints, solver!, constructor_beams, constructor_nodes, constructor_material_properties, constructor_geometry_properties, constructor_simulation_parameters, constructor_ext_force, constructor_boundary_conditions, constructor_configuration, Vec3, Vec2, Mat33, Geometry, Material, SDF_Plane_z, SDF_Plane_y, SDF_Sphere, SDF_Cylinder, constructor_discrete_sdf, read_TXT_file_pos, read_TXT_file_conn, read_TXT_file_ICs_array, read_TXT_file_ICs_matrix

end
