"""
Run a deployment simulation for a stent.

This function reconstructs the crimped stent from previous results, aligns it 
with the deployment origin, incorporates positioning displacements and rotations, 
and sets up the finite element model with constraints. It then runs the 
time-integration to simulate the deployment of the stent.

# Arguments
- `free_positions::Matrix{Float64}`: Initial nodal coordinates of the stent.
- `connectivity::Matrix{Int}`: Beam connectivity matrix.
- `constraints_connectivity::Matrix{Int}`: Connectivity matrix for penalty-based constraints.
- `deployment_origin_point::Vec3`: Reference 3D point for stent deployment (first node origin).
- `output_dir_crimping::String`: Directory containing crimping results.
- `output_dir_positioning::String`: Directory containing positioning results.
- `output_dir_deployment::String`: Directory to store deployment simulation results.

# Returns
The function exports the stent solution and  associated beam data to `output_dir_deployment`.
"""
function deployment(free_positions, connectivity, constraints_connectivity, deployment_origin_point, output_dir_crimping, output_dir_positioning, output_dir_deployment)
    
    #-------------------------------
    # Crimped stent
    #-------------------------------
    crimped_positions = matrix_to_vec3_array(free_positions) +
                        matrix_to_vec3_array(readdlm(output_dir_crimping * "u.txt"))

    # Compute initial stent centerline
    stent_centerline = compute_stent_centerline(crimped_positions)

    # Align stent first node with deployment origin
    shift_positions!(crimped_positions, stent_centerline[1] - deployment_origin_point)
    crimped_positions = vec3_array_to_matrix(crimped_positions)

    #-------------------------------
    # Positioned stent
    #-------------------------------
    positioned_positions = crimped_positions + readdlm(output_dir_positioning * "u.txt")
    initial_displacements = positioned_positions - free_positions 

    #-------------------------------
    #  Nodes
    #-------------------------------
    num_nodes = size(free_positions, 1)

    # Set initial rotation matrices for nodes
    initial_R⁰ = readdlm(output_dir_positioning * "R.txt")    

    # Create nodes structure    
    zeros_6d = zeros(size(free_positions))
    nodes = NodesBeams(
        free_positions,
        initial_displacements, zeros_6d, zeros_6d,  # zero initial displacements, velocities, accelerations
        zeros_6d, zeros_6d, zeros_6d,  # zero initial  rotations, angular velocities and accelerations
        nothing, initial_R⁰)

    #-------------------------------
    # Beams and constraints
    #-------------------------------
    E = 225*1e3
    ν = 0.33
    mass_scaling = 1e3
    ρ = 9.13*1e-9 * mass_scaling
    radius = 0.014
    damping = 1e4

    # Set initial rotation matrices for beams
    initial_Re⁰ = readdlm(output_dir_positioning * "Re0.txt")

    # Create beams structure
    beams = Beams(nodes, connectivity, E, ν, ρ, radius, damping, initial_Re⁰)
    num_beams = length(beams)

    # Create penalty constraints
    k_penalty = 1e3
    damping_penalty = 1
    constraints = Constraints(constraints_connectivity, k_penalty, damping_penalty)

    #-------------------------------
    # Assembled configuration
    #-------------------------------
    conf = BeamsConfiguration(nodes, beams, constraints, nothing, nothing)
    
    #-------------------------------
    # Simulation parameters
    #-------------------------------
    # HHT (Hilber–Hughes–Taylor) α-method parameters
    α = -0.05
    β = 0.25 * (1 - α)^2
    γ = 0.5 * (1 - 2 * α)
    
    params = SimulationParams(
        α = α, β = β, γ = γ,
        initial_timestep = 1e-6,
        min_timestep = 1e-9,
        max_timestep = 1e-3,
        output_timestep = 1e-3,
        simulation_end_time = 1,
        tolerance_residual = 1e-5,
        tolerance_displacement = 1e-5,
        max_iterations = 10,
        output_dir = output_dir_deployment,
        verbose = true
    )

    #-------------------------------
    # Run simulation and export results
    #-------------------------------
    run_simulation!(conf, params)
    export_stent_solution_to_txt(nodes, beams, num_nodes, num_beams, output_dir_deployment, false)

end 
