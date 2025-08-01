function positioning(initial_positions_stent, connectivity_stent, constraints_connectivity_stent, nb_iterations, origin, output_dir_crimping, output_dir_positioning_cl, output_dir_positioning)
    
    # --------
    # Stent 
    # --------
    crimped_positions = matrix_to_vec3_array(initial_positions_stent) +
                        matrix_to_vec3_array(readdlm(output_dir_crimping * "u.txt"))

    positions_cl =  compute_stent_centerline(crimped_positions)
    shift_positions!(crimped_positions, positions_cl[1]-origin)
    positions_cl =  compute_stent_centerline(crimped_positions, origin)

    # --------
    # Guides 
    # --------
    
    positions_guides, connectivity_guides = build_int_guides_stent_origin(crimped_positions, positions_cl[1])

    # ----------
    # Guides = stent
    # ------------

    initial_positions_stent = [crimped_positions; positions_guides]
    connectivity = [matrix_to_vec2_array(connectivity_stent); connectivity_guides]

    nnodes_stent = length(crimped_positions)
    nbeams_stent = length(matrix_to_vec2_array(connectivity_stent))
    nnodes = length(initial_positions_stent)  
    nbeams = length(connectivity)  
    iguides = nnodes_stent+1:nnodes

    # ----------
    # Building the nodes 
    # ----------
    
    # total number of nodes
    nnodes = size(initial_positions_stent, 1)
    
    # initial conditions
    initial_displacements = zeros(Vec3, nnodes)
    initial_velocities = zeros(Vec3, nnodes)
    initial_accelerations = zeros(Vec3, nnodes)
    initial_rotations = zeros(Vec3, nnodes)
    initial_angular_velocities = zeros(Vec3, nnodes)
    initial_angular_accelerations = zeros(Vec3, nnodes)
    plane = "xy"

    # Read nodes initial rotations
    R⁰ =  zeros(Mat33, nnodes)
    R⁰[1:nnodes_stent] .= matrix_to_mat33_array(readdlm(output_dir_crimping * "R.txt"))
    R⁰[nnodes_stent+1:nnodes] .= (Diagonal(Vec3(1,1,1)),)
    
    # Build nodes
    nodes = NodesBeams(initial_positions_stent, initial_displacements, initial_velocities, initial_accelerations, initial_rotations, initial_angular_velocities, initial_angular_accelerations, plane, R⁰)
    
    # -------
    # Building the beams
    # --------
    
    # Geometric and material properties
    E = 225*1e3
    ν = 0.33
    ρ = 9.13*1e-3
    radius = 0.065
    damping = 1E2

    # Read beams initial rotations
    Re₀ =  zeros(Mat33, nbeams)
    Re₀[1:nbeams_stent] .= matrix_to_mat33_array(readdlm(output_dir_crimping * "Re0.txt"))
    for (i,c) in enumerate(connectivity_guides)
        node1 = nodes[c[1]]
        node2 = nodes[c[2]]
        Re₀[nbeams_stent+i] = get_Rₑ⁰(node1.X₀, node2.X₀)
    end 

    # beams vector
    beams = Beams(nodes, connectivity, E, ν, ρ, radius, damping, Re₀)
    nbeams = length(beams)

    #----------------------------------
    # BEAMS CONFIGURATION DEFINITIONS
    #----------------------------------
    
    # Initialize the loads as nothing
    loads = nothing

    # number of dof (6 per node)
    ndofs = nnodes * 6
    
    # penalty constraints
    kᶜᵒⁿ = 1e3
    ηᶜᵒⁿ = 1
    constraints = Constraints(constraints_connectivity_stent, kᶜᵒⁿ, ηᶜᵒⁿ)
    
    # # Dirichlet boundary conditions: blocked positions
    # fixed_dofs = Vector{Int}() 
    # for i in iguides
    #     push!(fixed_dofs, 6*(i-1)+1)
    #     push!(fixed_dofs, 6*(i-1)+2)
    #     push!(fixed_dofs, 6*(i-1)+3)
    # end
    # free_dofs = setdiff(1:ndofs, fixed_dofs)
    
    # # Dirichlet dof (x6)
    # flag_load_from_file = true
    # dir_folder_load = output_dir_positioning_cl
    # disp_dofs = fixed_dofs
    # disp_vals = zeros(ndofs)
    # disp(t, node_idx) =  0

    # # boundary conditions strucutre
    # bcs = BoundaryConditions(fixed_dofs, free_dofs, disp, disp_vals, disp_dofs, flag_cylindrical, flag_load_from_file, dir_folder_load)
    bcs = nothing 

    # --------
    # Final configuration
    # ----------
    
    # configuration: mesh, external forces and boundary conditions
    conf = BeamsConfiguration(nodes, beams, constraints, loads, bcs)
    
    # ---------
    # Time stepping parameters
    # ---------
    
    # HHT (Houbolt-Hughes-Taylor) time stepping parameters
    α = -0.05 # Typically between 0 and -1, used for numerical stability
    β = 0.25 * (1 - α)^2 # Damping parameter based on α
    γ = 0.5 * (1 - 2 * α) # Time-stepping parameter
    
    # General time stepping parameters
    initial_timestep = 1 # Initial time step size
    min_timestep = 1 # Minimum allowed time step
    max_timestep = 1 # Maximum allowed time step (could be adjusted based on system behavior)
    output_timestep = 1 # Time step for output plotting or visualization
    simulation_end_time = nb_iterations-1 # End time for the simulation (duration of the analysis)
    
    # Convergence criteria for the solver
    tolerance_residual = 1e-5 # Residual tolerance for convergence checks
    tolerance_displacement = 1e-5 # Tolerance for changes in displacement (ΔD)
    max_iterations = 10 # Maximum number of iterations for the solver
    
    # Store solver parameters in a structured Params object
    params = SimulationParams(;
    α, β, γ, initial_timestep, min_timestep, max_timestep, output_timestep, simulation_end_time, tolerance_residual, tolerance_displacement, max_iterations, output_dir = output_dir_positioning, verbose = true)

    # -------------------
    # Start simulation
    # -------------------
    
    run_simulation!(conf, params)
    export_stent_solution_to_txt(nodes, beams, nnodes, nbeams, output_dir_crimping, false)

end 
