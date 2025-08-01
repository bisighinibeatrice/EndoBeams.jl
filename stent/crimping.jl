function crimping(rStent, positions, connectivity, constraints_connectivity, output_dir_crimping)
    
    #----------------------------------
    # BUILDING NODES 
    #----------------------------------
    
    # Number of nodes
    nnodes = size(positions, 1)
    
    # Initial conditions for displacements, velocities, accelerations, and rotations
    initial_displacements = zeros(size(positions)) # Zero initial displacements
    initial_velocities = zeros(size(positions)) # Zero initial velocities
    initial_accelerations = zeros(size(positions)) # Zero initial accelerations
    initial_rotations = zeros(size(positions)) # Zero initial rotations
    initial_angular_velocities = zeros(size(positions)) # Zero initial angular velocities
    initial_angular_accelerations = zeros(size(positions)) # Zero initial angular accelerations
    plane = "xy" # Plane of the problem
    
    # Build nodes
    nodes = NodesBeams(positions, initial_displacements, initial_velocities, initial_accelerations, initial_rotations, initial_angular_velocities, initial_angular_accelerations, plane) 
    
    #----------------------------------
    # BUILDING BEAMS AND CONSTRAINTS
    #----------------------------------
    
    # Geometric and material properties
    E = 225*1e3
    ν = 0.33
    mass_scaling = 1E3
    ρ = 9.13*1e-9 * mass_scaling
    radius = 0.014
    damping = 1E4
    
    # Build beams
    beams = Beams(nodes, connectivity, E, ν, ρ, radius, damping)
    nbeams = length(beams)
    
    # penalty constraints
    kᶜᵒⁿ = 1e3
    ηᶜᵒⁿ = 1
    constraints = Constraints(constraints_connectivity, kᶜᵒⁿ, ηᶜᵒⁿ)
    
    #----------------------------------
    # BEAMS CONFIGURATION DEFINITIONS
    #----------------------------------
    
    # Initialize the loads as nothing
    loads = nothing

    # Encastre
    ndofs = nnodes * 6                     # Total number of DOFs (6 per node for displacement and rotation)
    blocked_dofs = 2:6:ndofs-4              
    encastre = Encastre(blocked_dofs)

    # Imposed diplacement
    flag_cylindrical = true
    displaced_dof = 1:6:ndofs-5 # DOFs that will have the imposed displacement
    max_displacement = -0.5*rStent # Maximum displacement
    time_threshold = 1 # Time threshold for displacement change
    velocity = max_displacement / time_threshold # Velocity coefficient
    
    displacement_function(t) = (velocity*t)*(t<time_threshold) + (velocity*time_threshold)*(t>=time_threshold)
    
    imposed_displacement = ImposedDisplacement(displaced_dof, displacement_function)
    
    # Beam configuration struct initialization
    bcs =  BoundaryConditions(encastre, imposed_displacement, ndofs, flag_cylindrical)

    # bcs = nothing
    
    # Mesh configuration initialization
    conf = BeamsConfiguration(nodes, beams, constraints, loads, bcs)
    
    #----------------------------------------------------
    # SOLVER DEFINITIONS
    #----------------------------------------------------
    
    # HHT (Houbolt-Hughes-Taylor) time stepping parameters
    α = -0.05 # Typically between 0 and -1, used for numerical stability
    β = 0.25 * (1 - α)^2 # Damping parameter based on α
    γ = 0.5 * (1 - 2 * α) # Time-stepping parameter
    
    # General time stepping parameters
    initial_timestep = 1e-1 # Initial time step size
    min_timestep = 1e-6 # Minimum allowed time step
    max_timestep = 1e-1 # Maximum allowed time step (could be adjusted based on system behavior)
    output_timestep = 1e-1 # Time step for output plotting or visualization
    simulation_end_time = 1 # End time for the simulation (duration of the analysis)
    
    # Convergence criteria for the solver
    tolerance_residual = 1e-5 # Residual tolerance for convergence checks
    tolerance_displacement = 1e-5 # Tolerance for changes in displacement (ΔD)
    max_iterations = 10 # Maximum number of iterations for the solver
    
    # Store solver parameters in a structured Params object
    params = SimulationParams(;
    α, β, γ, initial_timestep, min_timestep, max_timestep, output_timestep, simulation_end_time, tolerance_residual, tolerance_displacement, max_iterations, output_dir = output_dir_crimping
    , verbose = true)
    
    #----------------------------------------------------
    # START SIMULATION
    #----------------------------------------------------
    
    # Run the solver with the defined configurations and parameters
    run_simulation!(conf, params)
    
end 

