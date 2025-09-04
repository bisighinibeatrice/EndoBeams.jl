# Function to create connectivity between stent ring nodes and the guide centerline
function get_guides_connectivity(rings, stent_centerline_indices)
    # Ensure each ring has a corresponding centerline node
    if length(rings) != length(stent_centerline_indices)
        error("length(rings) != length(stent_centerline).jl")
    end 

    nrings = length(rings)  # Total number of rings
    guides = Vector{Vec2{Int}}()  # Stores connections between ring nodes and centerline nodes
    
    # Loop over each ring
    for i in 1:nrings
        ring = rings[i]
        
        # Connect each node in the ring to the corresponding centerline node
        for node in ring
            push!(guides, Vec2(node, stent_centerline_indices[i]))
        end 
    end     
    
    return guides
end


# Function to build stent guides starting from an origin
function build_guides_stent_origin(positions_stent, origin)
    nnodes_stent = length(positions_stent)  # Total number of stent nodes
    
    # Group stent nodes by their z-coordinate to form rings
    rings = group_nodes_by_z_level(positions_stent)

    # Create the guide centerline nodes
    stent_centerline = Vector{Vec3{Float64}}()
    for ring in rings
        # Each centerline node inherits z from the first node of the ring
        push!(stent_centerline, [origin[1], origin[2], positions_stent[ring[1]][3]])
    end

    nnodes_stent_centerline = length(stent_centerline)

    # Assign indices to centerline nodes starting after the stent nodes
    stent_centerline_indices = nnodes_stent+1 : nnodes_stent + nnodes_stent_centerline

    # Build connectivity (these connections are called guides)
    connectivity_guides = get_guides_connectivity(rings, stent_centerline_indices)

    return stent_centerline, connectivity_guides
end
