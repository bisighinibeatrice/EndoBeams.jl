"""
Run a crimping simulation on a stent.

This function sets up the finite element model, boundary conditions, material 
properties, and time-integration scheme for simulating the crimping process 
of a stent represented as a network of beam elements.

# Arguments
- `free_radius::Float64`: Radius of the stent.
- `free_positions::Matrix{Float64}`: Initial nodal coordinates matrix.
- `connectivity::Matrix{Int}`: Beam connectivity matrix.
- `constraints_connectivity::Matrix{Int}`: Connectivity matrix for penalty-based constraints.
- `output_dir_crimping::String`: Directory to store simulation results.

# Returns
The function exports the stent solution and  associated beam data to `output_dir_crimping`.
"""
function crimping(free_radius, free_positions, connectivity, constraints_connectivity, output_dir_crimping)
    
    #-------------------------------
    # Nodes
    #-------------------------------
    num_nodes = size(free_positions, 1)

    # Create nodes structure
    zeros_6d = zeros(size(free_positions))  
    nodes = NodesBeams(
        free_positions,
        zeros_6d, zeros_6d, zeros_6d,  # zero initial displacements, velocities, accelerations
        zeros_6d, zeros_6d, zeros_6d,  # zero initial rotations, angular velocities, accelerations
        "xy"                           # working plane for the model (for cylindrical bcs)
    )
    
    #-------------------------------
    # Beams and constraints
    #-------------------------------
    E = 225e3                      # Young’s modulus (MPa) 
    ν = 0.33                       # Poisson’s ratio
    radius = 0.014                 # Beam radius (mm)
    ρ = 9.13e-9 * 1e3              # Density  (Tonne/mm3)
    damping = 1e4                  # Rayleigh damping coefficient
    
    # Create beams structure
    beams = Beams(nodes, connectivity, E, ν, ρ, radius, damping)
    num_beams = length(beams)

    # Create penalty constraints
    k_penalty = 1e3
    damping_penalty = 1
    constraints = Constraints(constraints_connectivity, k_penalty, damping_penalty)
    
    #-------------------------------
    # Loads and boundary conditions
    #-------------------------------
    num_dofs = num_nodes * 6

    # Fixed boundary condition (encastre)
    blocked_dofs = 2:6:num_dofs-4              
    encastre = Encastre(blocked_dofs)

    # Imposed displacement condition (radial crimping)
    max_disp = -0.7 * free_radius           # Inward radial displacement
    t_thresh = 1.0                          # Duration over which displacement is applied
    velocity = max_disp / t_thresh          # Constant displacement rate
    displaced_dof = 1:6:num_dofs-5
    displacement_fn(t) = velocity * min(t, t_thresh)

    imposed_disp = ImposedDisplacement(displaced_dof, displacement_fn)

    # Boundary conditions container
    bcs = BoundaryConditions(encastre, imposed_disp, num_dofs, true)

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
        initial_timestep = 1e-1,
        min_timestep = 1e-6,
        max_timestep = 1e-1,
        output_timestep = 1e-1,
        simulation_end_time = 1.0,
        tolerance_residual = 1e-5,
        tolerance_displacement = 1e-5,
        max_iterations = 10,
        output_dir = output_dir_crimping,
        verbose = true
    )

    #-------------------------------
    # Run simulation and export results
    #-------------------------------
    run_simulation!(conf, params)
    
    export_stent_solution_to_txt(
        nodes, beams, num_nodes, num_beams,
        output_dir_crimping,
        true
    )
end
