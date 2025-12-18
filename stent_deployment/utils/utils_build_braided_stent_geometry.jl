# Struct to define the parameters of a braided stent
@with_kw struct BraidedStent
    nbWires::Int = 12                      # Total number of wires used in the stent
    rStent::Float64 = 2.6                  # Radius of the deployed stent
    rCrimpedStent::Float64 = 1.3           # Radius when the stent is crimped
    rWireSection::Float64 = 0.051          # Radius of the individual wire cross-section
    wireGap::Float64 = 0                   # Optional gap between wires (zero = touching)
    lengthStent::Float64 = 5            # Length of the stent
    nbTotalCells::Float64 = 10             # Number of braid cells along the stent length
    braidingPattern::Int = 2               # Controls alternation pattern in braiding
end

"""
Generate a helical wire path (1 wire) and update positions and connectivity.

Arguments:
- `positions`: Vector of node positions (Vec3)
- `connectivity`: Wire connectivity (Vec2, node pairs)
- `nodeID`: Keeps track of node indices
- `length`, `diameter`: Stent dimensions
- `nbTotalCells`, `nbWires`: Braid resolution and density
- `rWireSection`, `braidingPattern`, `wireGap`: Geometry details
- `wireID`: Index of current wire
- `clockwise`: Boolean flag for wire direction
"""
function generate_helix_wire!(
    positions, connectivity, nodeID,
    length, diameter, nbTotalCells,
    nbWires, rWireSection, braidingPattern,
    wireGap, wireID, clockwise
)
    pitch = length / nbTotalCells

    orient = clockwise ? 1 : -1
    alternateBraiding = clockwise ? 0 : 1

    dTheta = 2π / (2 * nbWires)                  # Angular step between nodes
    offsetTheta = 2π / nbWires * wireID          # Angular offset for this wire

    braidingAngle = rad2deg(atan(pitch / (dTheta * diameter))) * 2
    nbZ = round(Int, 2 * nbTotalCells + 1)       # Total number of nodes along Z

    offsetBraiding = wireID % 2 == alternateBraiding ? -1 : 1

    for i in 0:(nbZ - 1)
        push!(nodeID, isempty(nodeID) ? 1 : nodeID[end] + 1)

        # Flip offset every 'braidingPattern' steps
        if (i % braidingPattern) == alternateBraiding
            offsetBraiding *= -1
        end

        # Calculate radial position
        r = diameter / 2
        if wireGap != 0
            r += offsetBraiding * orient * (rWireSection + wireGap / 2)
        end

        θ = orient * i * dTheta + offsetTheta
        x = r * cos(θ)
        y = r * sin(θ)
        z = i * pitch / 2

        push!(positions, [x, y, z])

        # Create edge (connect to previous node)
        if i != 0
            push!(connectivity, [nodeID[end - 1], nodeID[end]])
        end
    end
end

"""
Construct full braided stent geometry.

Returns:
- positions: Vector of 3D node positions
- connectivity: Vector of node pairs (wire segments)
"""
function generate_braided_stent_geometry(
    nbWires, radius, rWireSection, wireGap,
    length, nbTotalCells, braidingPattern
)
    positions = Vector{Vec3{Float64}}()
    connectivity = Vector{Vec2{Int}}()
    nodeID = Vector{Int}()

    # Create clockwise wires
    for wire in 1:nbWires
        generate_helix_wire!(positions, connectivity, nodeID, length, radius * 2,
               nbTotalCells, nbWires, rWireSection, braidingPattern,
               wireGap, wire, true)
    end

    # Create counter-clockwise wires
    for wire in 1:nbWires
        generate_helix_wire!(positions, connectivity, nodeID, length, radius * 2,
               nbTotalCells, nbWires, rWireSection, braidingPattern,
               wireGap, wire, false)
    end

    return positions, connectivity
end

"""
Group nodes by their Z-coordinate into rings.

Returns:
- `rings`: Vector of rings, each ring is a vector of node indices
"""
function group_nodes_by_z_level(positions)
    zvec = [pos[3] for pos in positions]         # Extract Z values
    zvec .= round.(zvec, digits = 4)             # Round to remove floating-point noise
    unique!(zvec)

    rings = Vector{Vector{Int}}()
    for z in zvec
        ring = Int[]
        for (i, p) in enumerate(positions)
            if isapprox(round(p[3], digits = 4), z)
                push!(ring, i)
            end
        end
        push!(rings, ring)
    end

    return rings
end

"""
Identify node pairs (within each ring) that are close enough to form constraints.

Returns:
- `constraints`: Vector of Vec2{Int} representing node pairs
"""
function find_close_node_pairs_in_rings(positions)
    constraints = Vector{Vec2{Int}}()
    tolerance = 0.06  # Distance threshold for creating constraint
    rings = group_nodes_by_z_level(positions)

    for ring in rings
        for i in 1:length(ring)
            idx_i = ring[i]
            pos_i = positions[idx_i]
            for j in 1:length(ring)
                idx_j = ring[j]
                pos_j = positions[idx_j]

                if idx_i != idx_j && norm(pos_i - pos_j) < tolerance
                    # Avoid duplicate pairs
                    if !([idx_j, idx_i] in constraints)
                        push!(constraints, [idx_i, idx_j])
                    end
                end
            end
        end
    end

    return constraints
end
