"""
Run crimping simulation on a stent using beam-based structural mechanics.

# Arguments
- `r_stent::Float64`: Target radius of the crimped stent.
- `initial_positions_stent::Matrix{Float64}`: Node initial_positions_stent.
- `connectivity_stent::Matrix{Int}`: Beam connectivity_stent (pairs of node indices).
- `constraints_connectivity_stent::Matrix{Int}`: connectivity_stent for penalty constraints.
- `output_dir_crimping::String`: Directory to store simulation results.
"""
function crimping(r_stent, initial_positions_stent, connectivity_stent, constraints_connectivity_stent, output_dir_crimping)
    
    #-------------------------------
    # NODE INITIALIZATION
    #-------------------------------
    
    nnodes = size(initial_positions_stent, 1)
    zeros_6d = zeros(size(initial_positions_stent))
    
    nodes = NodesBeams(
        initial_positions_stent,
        zeros_6d, zeros_6d, zeros_6d,  # zero initial displacements, velocities, accelerations
        zeros_6d, zeros_6d, zeros_6d,  # zero initial  rotations, angular velocities and accelerations
        "xy"                           # working plane
    )
    
    #-------------------------------
    # BEAMS AND CONSTRAINTS
    #-------------------------------

    E = 225e3                      # Young's modulus (MPa)
    ν = 0.33                       # Poisson's ratio
    radius = 0.014                 # Radius of the beam (mm)
    ρ = 9.13e-9 * 1e3              # Density (scaled)
    damping = 1e4                  # Rayleigh damping
    beams = Beams(nodes, connectivity_stent, E, ν, ρ, radius, damping)

    constraints = Constraints(constraints_connectivity_stent, 1e3, 1.0)  # penalty: stiffness and damping

    #-------------------------------
    # BOUNDARY CONDITIONS
    #-------------------------------
    
    # Encastre
    ndofs = nnodes * 6                
    blocked_dofs = 2:6:ndofs-4              
    encastre = Encastre(blocked_dofs)

    # Imposed displacement setup (crimping)
    max_disp = -0.5 * r_stent
    t_thresh = 1.0
    velocity = max_disp / t_thresh
    displaced_dof = 1:6:ndofs-5
    displacement_fn(t) = velocity * min(t, t_thresh)

    imposed_disp = ImposedDisplacement(displaced_dof, displacement_fn)

    # Beam configuration struct initialization
    bcs = BoundaryConditions(encastre, imposed_disp, ndofs, true)

    # Mesh configuration initialization
    conf = BeamsConfiguration(nodes, beams, constraints, nothing, bcs)
    
    #-------------------------------
    # CONFIGURATION 
    #-------------------------------

    conf = BeamsConfiguration(nodes, beams, constraints, nothing, bcs)

    #-------------------------------
    # PARAMETERS 
    #-------------------------------

    # HHT (Houbolt-Hughes-Taylor) time stepping parameters
    α = -0.05 # Typically between 0 and -1, used for numerical stability
    β = 0.25 * (1 - α)^2 # Damping parameter based on α
    γ = 0.5 * (1 - 2 * α) # Time-stepping parameter
    
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
    # RUN SIMULATION
    #-------------------------------

    run_simulation!(conf, params)
    export_stent_solution_to_txt(nodes, beams, nnodes, length(beams), output_dir_crimping, false)

end 

