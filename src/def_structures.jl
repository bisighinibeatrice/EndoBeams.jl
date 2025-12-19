# Structure for force vectors required at each time step in the solver
struct Forces
    fᵉˣᵗ::Vector{Float64}    # External forces vector
    Tⁱⁿᵗ::Vector{Float64}    # Internal forces vector
    Tᵏ::Vector{Float64}      # Kinetic forces vector (inertia effects)
    Tᶜ::Vector{Float64}      # Contact forces vector
end  

# Constructor for Solution structure, initializing force vectors based on beams configuration
function Forces(conf::BeamsConfiguration)

    ndofs = conf.ndofs  # Total dofs

    # Initialize external forces vector with zeros
    fᵉˣᵗ = zeros(ndofs)

    # Populate external forces vector based on loaded dofs in configuration
    if conf.loads !== nothing && conf.loads.concentrated_force !== nothing
        for i in conf.loads.concentrated_force.loaded_dofs
            fᵉˣᵗ[i] = conf.loads.concentrated_force.force_function(0,i)
        end
    end 
    
    # Initialize remaining force vectors with zeros
    Tⁱⁿᵗ = zeros(ndofs)  # Internal forces
    Tᵏ = zeros(ndofs)     # Kinetic forces
    Tᶜ = zeros(ndofs)     # Contact forces

    return Forces(fᵉˣᵗ, Tⁱⁿᵗ, Tᵏ, Tᶜ)
end

# Structure to hold nodal solutions with preallocated vectors for solver computations
struct Solution 
    D::Vector{Float64}              # Displacement vector
    Ḋ::Vector{Float64}             # Velocity vector
    D̈::Vector{Float64}             # Acceleration vector
    Ḋⁿ::Vector{Float64}             # Old velocity vector
    r::Vector{Float64}              # Residual forces vector
    Ktan::SparseMatrixCSC{Float64,Int} # Global stiffness matrix in sparse format
    ΔD::Vector{Float64}             # Displacement increment vector
    temp::Vector{Float64}           # Temporary vector for intermediate computations
    r_free::Vector{Float64}         # Free dofs residual vector
    Ktan_free::SparseMatrixCSC{Float64,Int} # Free DOF stiffness matrix
    ΔD_free::Vector{Float64}        # Free DOF displacement increment vector
end 

# Constructor for Solution structure, setting up required vector and matrix allocations
function Solution(Ktan, Ktan_free, ndofs, nfreedofs)
    return Solution(
        zeros(ndofs),          # Displacement vector
        zeros(ndofs),          # Velocity vector
        zeros(ndofs),          # Acceleration vector
        zeros(ndofs),          # Old velocity vector
        zeros(ndofs),          # Residual forces vector
        Ktan,                  # Global stiffness matrix (sparse)
        zeros(ndofs),          # Displacement increment
        zeros(ndofs),          # Temporary vector
        zeros(nfreedofs),      # Free DOF residual forces vector
        Ktan_free,             # Free DOF stiffness matrix
        zeros(nfreedofs)       # Free DOF displacement increment
    )
end

# Structure to hold global matrices (stiffness, damping, mass) for solver operations
struct Matrices
    K::SparseMatrixCSC{Float64,Int} # Global stiffness matrix (sparse)
    C::SparseMatrixCSC{Float64,Int} # Damping matrix (sparse)
    M::SparseMatrixCSC{Float64,Int} # Mass matrix (sparse)
    sparsity_free::Vector{Int}          # Free dofs sparsity pattern
end 

# Constructor for Matrices structure, initializing matrices and forces for solver computations
function Matrices(I, J, sparsity_free)

    # Create sparse matrices for stiffness, damping, and mass
    K = sparse(I, J, 0.)
    C = sparse(I, J, 0.)
    M = sparse(I, J, 0.)
    
    return Matrices(K, C, M, sparsity_free)
end

# Structure to store energy contributions in the simulation
mutable struct Energy
    strain_energy::Float64      # Strain energy contribution
    kinetic_energy::Float64     # Kinetic energy contribution
    contact_energy::Float64     # Contact energy contribution
end

# Constructor for Energy structure, initializing energy values to zero
function Energy()
    return Energy(0.0, 0.0, 0.0)
end

# Structure to store the state of the system at different time steps
struct SimulationState 
    forcesⁿ::Forces          # Forces at the current time step for the system
    forcesⁿ⁺¹::Forces        # Forces at the next time step  for the system
    matricesⁿ::Matrices      # Matrices at the next time step  for the system
    matricesⁿ⁺¹::Matrices    # Matrices at the current time step for the system
    solⁿ⁺¹::Solution         # Solution for the system
    energyⁿ⁺¹::Energy        # Energy state for the system
end 
