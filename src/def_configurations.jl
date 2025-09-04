# Abstract type for mesh configurations
abstract type SimulationConfiguration end

# Structure defining a beam mesh configuration and its properties
struct BeamsConfiguration{Tn, Tb, Tcon} <: SimulationConfiguration
    nodes::Tn                 # Mesh nodes
    beams::Tb                 # Beam elements
    ndofs::Int                # Total degrees of freedom (DOFs)
    offset_dofs::Int          # Offset DOFs (always 0 since beams are allocated first if solids exist -not in this version of the code)
    constraints::Tcon         # Constraints
    loads::Union{Loads, Nothing}           # External forces applied to the mesh
    bcs::Union{BoundaryConditions, Nothing}                  # Boundary conditions
end

# Constructor for BeamsConfiguration (for beam meshes)
function BeamsConfiguration(nodes::StructVector, beams::StructVector, constraints::Union{StructVector, Nothing}, loads::Union{Loads, Nothing}, bcs::Union{BoundaryConditions, Nothing})
    
    ndofs = length(nodes) * 6  # Each node has 6 DOFs (3 displacements + 3 rotations)
    
    if bcs === nothing bcs = BoundaryConditions(ndofs) end 

    return BeamsConfiguration{typeof(nodes), typeof(beams), typeof(constraints)}(
    nodes, beams, ndofs, 0, constraints, loads, bcs)

end

# Constructor for BeamsConfiguration (for beam meshes) without constraints
function BeamsConfiguration(nodes::StructVector, beams::StructVector, loads::Union{Loads, Nothing}, bcs::Union{BoundaryConditions, Nothing})
    
    ndofs = length(nodes) * 6  # Each node has 6 DOFs (3 displacements + 3 rotations)
    
    if bcs === nothing bcs = BoundaryConditions(ndofs) end 

    constraints = nothing 
    
    return BeamsConfiguration{typeof(nodes), typeof(beams), typeof(constraints)}(
    nodes, beams, ndofs, 0, nothing, loads, bcs)

end