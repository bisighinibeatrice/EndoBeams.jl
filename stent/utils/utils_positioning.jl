
function get_internal_nodes(positions, rStent, rWireSection, wireGap)
    
    toll = (rWireSection + wireGap/2)/2
    
    Rint = rStent - rWireSection - wireGap/2
    Rext = rStent + rWireSection + wireGap/2
    
    intnodes = []
    extnodes = []
    
    for i in 1:length(positions)
        
        r = sqrt(positions[i][1]^2 + positions[i][2]^2)
        
        if isapprox(r, Rint; atol = toll) 
            push!(intnodes, i)
        elseif isapprox(r, Rext; atol = toll) 
            push!(extnodes, i)
        end 
        
    end 
    
    return intnodes, extnodes
    
end

function get_guide_connectivity_int(rings, nodes_centerline)
    
    if length(rings) != length(nodes_centerline)
        error("length(rings) != length(nodes_centerline).jl")
    end 

    nrings = length(rings)
    constraints = Vector{Vec2{Int}}()
    
    for i in 1:nrings
        thisRing = rings[i]
        nnodesRing = length(thisRing)
        for j in 1:nnodesRing
            push!(constraints, Vec2(thisRing[j], nodes_centerline[i]))
        end 
    end     
    
    return constraints
    
end

function build_int_guides_stent_origin(positions, origin, type="Surpass")

    thisWireGap = 0.1

    # total number of stent nodes
    nnodes_stent = length(positions)
    
    # get rings
    rings = group_nodes_by_z_level(positions)

    #create centerline 
    positions_guides = Vector{Vec3{Float64}}()
    nguides = length(rings)

    for ring in 1:nguides
        push!(positions_guides, [origin[1], origin[2], positions[rings[ring][1]][3]])
    end
    nnodes_centerline = length(positions_guides)
    nodes_centerline = nnodes_stent+1:nnodes_stent + nnodes_centerline

    # # get internal node
    # if type == "Surpass"    
    #     @unpack nbWires, rStent, rWireSection, wireGap, lengthStent, nbTotalCells, braidingPattern = BraidedStent()
    #     positions, connectivity_stent = compute_bs_geom_given_nbTotalCells(nbWires, rStent, rWireSection, thisWireGap, lengthStent, nbTotalCells, braidingPattern)
    #     intnodes, extnodes = get_internal_nodes(positions, rStent, rWireSection, thisWireGap*1.1)
    # else type == "Wallstent" 
    #     @unpack nbWires, rStent, rCrimpedStent, rWireSection, lStent, phi, braidingPattern = Wallstent()
    #     positions, connectivity_stent = compute_bs_geom(nbWires, rStent, rCrimpedStent, rWireSection, thisWireGap, lStent, phi, braidingPattern)
    #     intnodes, extnodes = get_internal_nodes(positions, rStent, rWireSection, thisWireGap*1.1)
    # end

    # if length(intnodes) != length(extnodes)
    #     println("length(intnodes) : ")
    #     println(length(intnodes))
    #     println("length(extnodes) : ")
    #     println(length(extnodes))
    #     error("Error in the computation of the internal nodes. ")
    # end 

    connectivity_guides = get_guide_connectivity_int(rings, nodes_centerline)

    return positions_guides, connectivity_guides

end