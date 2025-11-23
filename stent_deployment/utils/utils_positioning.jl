"""
Create connectivity between stent ring nodes and the guide centerline.

Each node in a ring is connected to the corresponding centerline node 
at the same z-level.

# Arguments
- `rings::Vector{Vector{Int}}`: List of rings, where each ring is a vector of 
 stent node indices belonging to the same z-level.
- `centerline_indices::UnitRange{Int}`: Indices of the guide centerline 
 nodes corresponding to each ring.

# Returns
- `Vector{Vec2{Int}}`: Connectivity list, where each element is a `Vec2` 
 representing a connection between a ring node and a centerline node.
"""
function get_guides_connectivity(rings, centerline_indices)
    if length(rings) != length(centerline_indices)
        error("length(rings) != length(centerline_indices)")
    end 
    
    nrings = length(rings)
    guides = Vector{Vec2{Int}}()
    
    for i in 1:nrings
        ring = rings[i]
        for node in ring
            push!(guides, Vec2(node, centerline_indices[i]))
        end 
    end 
    
    return guides
end

"""
Build guide structures for a stent starting from a given deployment_origin_point.

The guides are formed by grouping stent nodes into rings at each z-level, 
creating new centerline nodes aligned with the stent deployment_origin_point, and connecting 
the stent rings to these centerline nodes.

# Arguments
- `positions::Vector{Vec3{Float64}}`: Positions of stent nodes in 3D.
- `deployment_origin_point::Vec3{Float64}`: Reference deployment_origin_point point for the guide centerline.

# Returns
- `centerline::Vector{Vec3{Float64}}`: List of new centerline nodes 
 created along the z-axis of the stent.
- `connectivity_guides::Vector{Vec2{Int}}`: Connectivity linking stent ring 
 nodes to their corresponding centerline nodes.
"""
function build_guides_stent_deployment_origin_point(positions, deployment_origin_point)
    num_nodes = length(positions)
    
    # Group stent nodes by z-level into rings
    rings = group_nodes_by_z_level(positions)
    
    # Create guide centerline nodes aligned with deployment_origin_point (x, y) but varying in z
    centerline = Vector{Vec3{Float64}}()
    for ring in rings
        push!(centerline, [deployment_origin_point[1], deployment_origin_point[2], positions[ring[1]][3]])
    end
    num_nodes_centerline = length(centerline)
    
    # Assign indices to new centerline nodes after the stent nodes
    centerline_indices = num_nodes+1 : num_nodes + num_nodes_centerline
    
    # Build connectivity linking stent rings to centerline nodes
    connectivity_guides = get_guides_connectivity(rings, centerline_indices)
    
    return centerline, connectivity_guides
end

"""
Shift all positions in `positions` by vector `disp`, in-place.

# Arguments
- `positions::Vector{Vec3}`: Positions to modify.
- `disp::Vec3`: Displacement to subtract (e.g., for recentering).

# Side Effects
- Modifies `positions` in-place.
"""
function shift_positions!(positions, disp::Vec3)
    for i in eachindex(positions)
        positions[i] -= disp
    end
end

"""
Compute a centerline along the stent by averaging the XY positions of each Z-level ring.

# Arguments
- `positions::Vector{Vec3}`: All node positions in the stent mesh.
- `origin::Vec3`: Optional fixed XY reference (defaults to (0,0,0)).

# Returns
- `centerline::Vector{Vec3}`: Centerline points, one per Z-ring.
"""
function compute_stent_centerline(positions, origin::Vec3 = Vec3(0, 0, 0))
    rings = group_nodes_by_z_level(positions)
    centerline = Vector{Vec3}()
    
    for ring in rings
        z_val = mean(p[3] for p in positions[ring])
        center_point = Vec3(origin[1], origin[2], z_val)
        push!(centerline, center_point)
    end
    
    return centerline
end