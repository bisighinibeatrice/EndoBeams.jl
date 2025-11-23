module EndoBeams

    include("packages.jl")
    include("constants.jl")
    include("linearsolvers.jl")
    include("def_params.jl") 

    include("constraints/def_constraint.jl")

    include("def_bcs_and_loads.jl")
    include("def_configurations.jl")
    include("def_structures.jl")
    include("interactions/def_interactions.jl")
    include("sparse_matrices.jl")

    include("beams/def_node_beam.jl")
    include("beams/def_beam.jl")
    include("beams/compute_beams.jl")
    include("beams/assemble_beams.jl")
    include("beams/utils_beams.jl")
    include("constraints/compute_and_assemble_constraints.jl")

    include("interactions/def_interactions.jl")
    include("interactions/compute&assemble_interactions_beams.jl")
    include("interactions/interactions_beams_regularization.jl")
    include("interactions/utils_compute_interactions_beams.jl")

    include("inizialization.jl")
    include("run_simulation.jl")
    include("solve_step_dynamics.jl")
    include("predictor.jl")
    include("corrector.jl")
    include("utils_bcs_and_loads.jl")
    include("utils_solver.jl")

    include("IO/visualization_beams.jl") 
    include("IO/visualization.jl") 
    include("IO/read_input_files.jl")  

    export NodesBeams, Beams, Constraints, BeamsConfiguration
    export ConcentratedForce, Loads
    export Encastre, ImposedDisplacement,ImposedDisplacementFromFile, BoundaryConditions
    export SimulationParams, run_simulation!
    export read_vtk_tetrahedral_mesh, read_vtk_triangle_mesh
    export RigidInteraction, SoftInteraction, PlaneSurface, SphereSurface, TriangulatedSurface, BeamElementSurface, BeamNodeSurface, InteractionProperties, DiscreteSignedDistanceField
    export Matrices
    export Vec2, Vec3, Mat33, get_Rₑ⁰
    
    include("precompile/precompiles.jl")
    _precompile_()

end 