"""
Run a positioning simulation for a stent along its deployed centerline.

This function reconstructs the crimped stent from previous simulation results, 
aligns it with the deployment origin, adds guide structures, sets up the finite 
element model including boundary conditions and constraints, and runs the 
time-integration to simulate the positioning of the stent along the target 
centerline.

# Arguments
- `free_positions::Matrix{Float64}`: Initial nodal coordinates of the stent.
- `connectivity::Matrix{Int}`: Connectivity matrix for the stent beams.
- `constraints_connectivity::Matrix{Int}`: Connectivity matrix for penalty-based constraints.
- `num_morphing_iterations::Int`: Number of morphing/positioning iterations.
- `deployment_origin_point::Vec3`: Reference 3D point for stent deployment (first node origin).
- `output_dir_crimping::String`: Directory containing crimping simulation results.
- `output_dir_centerline_morphing::String`: Directory containing morphing results (used for imposed displacement on guides).
- `output_dir_positioning::String`: Directory to store the positioning simulation results.

# Returns
The function exports the stent solution and  associated beam data to `output_dir_positioning`.
"""
function positioning(free_positions, connectivity, constraints_connectivity, num_morphing_iterations, deployment_origin_point, output_dir_crimping, output_dir_centerline_morphing, output_dir_positioning)

    #-------------------------------
    # Crimped stent
    #-------------------------------
    crimped_positions = matrix_to_vec3_array(free_positions) +
                        matrix_to_vec3_array(readdlm(output_dir_crimping * "u.txt"))

    # Compute initial stent centerline
    stent_centerline = compute_stent_centerline(crimped_positions)

    # Align stent first node with deployment origin
    shift_positions!(crimped_positions, stent_centerline[1] - deployment_origin_point)

    # Recompute centerline after alignment
    stent_centerline = compute_stent_centerline(crimped_positions, deployment_origin_point)

    #-------------------------------
    # Guides
    #-------------------------------
    guide_positions, guide_connectivity = build_guides_stent_origin(crimped_positions, stent_centerline[1])

    #-------------------------------
    # Combining stent and guides 
    #-------------------------------
    all_positions = [crimped_positions; guide_positions]
    all_connectivity = [matrix_to_vec2_array(connectivity); guide_connectivity]

    num_stent_nodes = length(crimped_positions)
    num_stent_beams = length(matrix_to_vec2_array(connectivity))
    num_total_nodes = length(all_positions)
    num_total_beams = length(all_connectivity)
    guide_node_indices = num_stent_nodes+1:num_total_nodes

    #-------------------------------
    #  Nodes
    #-------------------------------
    initial_displacements = zeros(Vec3, num_total_nodes)
    initial_velocities = zeros(Vec3, num_total_nodes)
    initial_accelerations = zeros(Vec3, num_total_nodes)
    initial_rotations = zeros(Vec3, num_total_nodes)
    initial_angular_velocities = zeros(Vec3, num_total_nodes)
    initial_angular_accelerations = zeros(Vec3, num_total_nodes)
    working_plane = "xy"

    # Set initial rotation matrices for stent and guide nodes
    initial_R⁰ = zeros(Mat33, num_total_nodes)
    initial_R⁰[1:num_stent_nodes] .= matrix_to_mat33_array(readdlm(output_dir_crimping * "R.txt"))
    initial_R⁰[num_stent_nodes+1:num_total_nodes] .= (Diagonal(Vec3(1,1,1)),)

    # Create nodes structure
    nodes = NodesBeams(
        all_positions,
        initial_displacements,
        initial_velocities,
        initial_accelerations,
        initial_rotations,
        initial_angular_velocities,
        initial_angular_accelerations,
        working_plane,
        initial_R⁰
    )

    #-------------------------------
    # Beams and constraints
    #-------------------------------
    E_modulus = 225e3
    poisson_ratio = 0.33
    mass_scaling = 1e6
    density = 9.13e-9 * mass_scaling
    beam_radius = 0.065
    damping_coefficient = 1e4

    # Set initial rotation matrices for beams
    initial_Re⁰ = zeros(Mat33, num_total_beams)
    initial_Re⁰[1:num_stent_beams] .= matrix_to_mat33_array(readdlm(output_dir_crimping * "Re0.txt"))
    for (i, conn) in enumerate(guide_connectivity)
        node1 = nodes[conn[1]]
        node2 = nodes[conn[2]]
        initial_Re⁰[num_stent_beams + i] = get_Rₑ⁰(node1.X₀, node2.X₀)
    end

    # Create beams structure
    beams = Beams(nodes, all_connectivity, E_modulus, poisson_ratio, density, beam_radius, damping_coefficient, initial_Re⁰)

    # Create penalty constraints
    k_penalty = 1e3
    damping_penalty = 1
    constraints = Constraints(constraints_connectivity, k_penalty, damping_penalty)

    #-------------------------------
    # Loads and boundary conditions
    #-------------------------------
    num_dofs = num_total_nodes * 6

    # Imposed displacement condition (read from file)
    displaced_dofs = vcat([6*(i-1) .+ (1:3) for i in guide_node_indices]...)
    imposed_displacement = ImposedDisplacementFromFile(displaced_dofs, output_dir_centerline_morphing)
    bcs = BoundaryConditions(imposed_displacement, num_dofs)
    
    #-------------------------------
    # Assembled configuration
    #-------------------------------
    conf = BeamsConfiguration(nodes, beams, constraints, nothing, bcs)

    #-------------------------------
    # Simulation parameters
    #-------------------------------
    # HHT (Hilber–Hughes–Taylor) α-method parameters
    α = -0.05
    β = 0.25 * (1 - α)^2
    γ = 0.5 * (1 - 2 * α)

    params = SimulationParams(
        α = α, β = β, γ = γ,
        initial_timestep = 1.0,
        min_timestep = 1e-5,
        max_timestep = 1.0,
        output_timestep = 1.0,
        simulation_end_time = num_morphing_iterations - 1,
        tolerance_residual = 1e-5,
        tolerance_displacement = 1e-5,
        max_iterations = 10,
        output_dir = output_dir_positioning,
        verbose = true
    )

    #-------------------------------
    # Run simulation and export results
    #-------------------------------
    run_simulation!(conf, params)

    export_stent_solution_to_txt(nodes, beams, num_stent_nodes, num_stent_beams, output_dir_positioning, false)
end
