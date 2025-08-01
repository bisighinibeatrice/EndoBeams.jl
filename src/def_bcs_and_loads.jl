#----------------------------------
# BUILDING BOUNDARY CONDITIONS 
#----------------------------------

# Structure representing encastre (fully fixed) boundary conditions
struct Encastre 
    blocked_dofs::Union{Int, StepRange{Int, Int}, UnitRange{Int}, Vector{Int}}  # Degrees of freedom (DOFs) that are fully fixed
end
 
# Structure representing imposed displacement boundary conditions
struct ImposedDisplacement{TD}
    displaced_dofs::Union{Int, StepRange{Int, Int}, UnitRange{Int}, Vector{Int}}  # DOFs where displacements are imposed
    imposed_displacements_function::TD  # Function defining imposed displacement values over time
end

# Structure representing general boundary conditions in the system
struct BoundaryConditions{TD} 
    fixed_dofs::Union{Int, StepRange{Int, Int}, UnitRange{Int}, Vector{Int}}  # DOFs that are fixed
    free_dofs::Union{Int, StepRange{Int, Int}, UnitRange{Int}, Vector{Int}}   # DOFs that are free to move
    displaced_indices::Union{Int, StepRange{Int, Int}, UnitRange{Int}, Vector{Int}}  # Indices of imposed displacement DOFs among the fixed DOFs
    imposed_displacements::Vector{Float64}  # Vector storing imposed displacement values
    imposed_displacements_function::TD  # Function defining imposed displacements
    flag_cylindrical::Bool # Flag to check if BCs are applied in carthesian or cylindrical coordinate system
end

# Constructor for BoundaryConditions with encastre conditions (fully fixed DOFs)
function BoundaryConditions(bcs::Encastre, ndofs, flag_cylindrical::Union{Bool, Nothing} = false) 
    fixed_dofs = bcs.blocked_dofs  # Get blocked DOFs
    free_dofs = setdiff(1:ndofs, fixed_dofs)  # Compute free DOFs
    ndofs_fixed = length(fixed_dofs)
    imposed_displacements = zeros(ndofs_fixed)  # No imposed displacements
    displaced_indices = Int[]  # No displaced indices for encastre conditions

    if isnothing(flag_cylindrical) flag_cylindrical = false end 

    return BoundaryConditions{Nothing}(fixed_dofs, free_dofs, displaced_indices, imposed_displacements, nothing, flag_cylindrical)
end

# Constructor for BoundaryConditions with imposed displacement conditions
function BoundaryConditions(bcs::ImposedDisplacement, ndofs, flag_cylindrical::Union{Bool, Nothing} = false) 
    fixed_dofs = bcs.displaced_dofs  # Get imposed displacement DOFs
    free_dofs = setdiff(1:ndofs, fixed_dofs)  # Compute free DOFs
    ndofs_fixed = length(fixed_dofs)
    imposed_displacements = zeros(ndofs_fixed)  # Initialize imposed displacement values
    displaced_indices = 1:length(fixed_dofs)  # Store indices of imposed displacement DOFs

    if isnothing(flag_cylindrical) flag_cylindrical = false end 

    return BoundaryConditions{typeof(bcs.imposed_displacements_function)}(fixed_dofs, free_dofs, displaced_indices, imposed_displacements, bcs.imposed_displacements_function, flag_cylindrical)
end 

# Constructor for BoundaryConditions combining encastre and imposed displacement conditions
function BoundaryConditions(bcs1::Encastre, bcs2::ImposedDisplacement, ndofs, flag_cylindrical::Union{Bool, Nothing} = false) 
    fixed_dofs = sort(vcat(bcs2.displaced_dofs, bcs1.blocked_dofs))  # Combine and sort all fixed DOFs
    free_dofs = setdiff(1:ndofs, fixed_dofs)  # Compute free DOFs
    ndofs_fixed = length(fixed_dofs)
    imposed_displacements = zeros(ndofs_fixed)  # Initialize imposed displacement values
    displaced_indices = findall(x -> x in bcs2.displaced_dofs, fixed_dofs)  # Find imposed displacement indices

    if isnothing(flag_cylindrical) flag_cylindrical = false end 

    return BoundaryConditions{typeof(bcs2.imposed_displacements_function)}(fixed_dofs, free_dofs, displaced_indices, imposed_displacements, bcs2.imposed_displacements_function, flag_cylindrical)
end

# Constructor for BoundaryConditions when no explicit boundary conditions are provided (all DOFs free)
function BoundaryConditions(ndofs) 
   
    fixed_dofs = Int[]  # No fixed DOFs (all free)
    free_dofs = 1:ndofs  # All DOFs are free
    imposed_displacements = Float64[]  # No imposed displacements
    displaced_indices = Int[]  # No imposed displacement indices

    return BoundaryConditions{Nothing}(fixed_dofs, free_dofs, displaced_indices, imposed_displacements, nothing, false)
end

#----------------------------------
# BUILDING LOADS 
#----------------------------------

# Structure representing a concentrated force applied at specific DOFs
struct ConcentratedForce{TF}
    force_function::TF  # Function defining applied force over time
    loaded_dofs::Vector{Int}  # DOFs where forces are applied
end

# Constructor for ConcentratedForce
function ConcentratedForce(force_fun::TF, loaded_dofs) where TF
    return ConcentratedForce{TF}(force_fun, loaded_dofs)
end

# Structure representing all applied loads in the model
struct Loads
    concentrated_force::ConcentratedForce  # Applied concentrated forces
    # distributed_force::DistributedForce  # To be added in the future
    # gravity::Gravity  # To be added in the future
end
