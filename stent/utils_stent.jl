# -------------------------------------------------------------------------------------------
# Stent
# -------------------------------------------------------------------------------------------

@with_kw struct BraidedStent 

    nbWires::Int = 24
    rStent::Float64 = 2.3
    rCrimpedStent::Float64 = 1.1
    rWireSection::Float64 = 0.014
    wireGap::Float64 = 0
    lengthStent::Float64 = 7.5
    nbTotalCells::Float64 = 35
    braidingPattern::Int = 2

end

function helix!(positions, connectivity, nodeID, length, diameter, nbTotalCells, nbWires, rWireSection, braidingPattern, wireGap, wireID, clockwise)
    
    pitch = length / nbTotalCells
    
    if clockwise
        orient = 1
        alternateBraiding = 0
    else
        orient = -1
        alternateBraiding = 1
    end 
    
    dTheta = 2*pi/(2*nbWires)
    
    offsetTheta = 2*pi/(nbWires)*wireID
    braidingAngle = 360/(2*pi)*tan(pitch/(dTheta*diameter))*2
    invBraidingAngle = 180 - braidingAngle
    nbZ = round(2*nbTotalCells+1)
    
    if wireID%2 == alternateBraiding
        offsetBraiding = -1
    else
        offsetBraiding = 1 
    end 
    
    for i in 0:nbZ-1
        
        if nodeID == []
            push!(nodeID, 1)
        else 
            push!(nodeID, nodeID[end]+1)
        end 
        
        if (i%braidingPattern) == alternateBraiding
            offsetBraiding = -offsetBraiding
        end 
        
        if wireGap != 0
            r = diameter/2 + offsetBraiding * orient * (rWireSection+wireGap/2)
        else 
            r = diameter/2
        end 
        x_n = r*cos(orient*i*dTheta+offsetTheta)
        y_n = r*sin(orient*i*dTheta+offsetTheta)
        z_n = i*pitch/2
        
        push!(positions, [x_n, y_n, z_n])
        
        if i != 0 
            push!(connectivity, [nodeID[end-1], nodeID[end]])
        end 
    end
 
end 

function compute_bs_geom_given_nbTotalCells(nbWires, radius, rWireSection, wireGap, length, nbTotalCells, braidingPattern)
    
    positions = Vector{Vec3{Float64}}()
    connectivity = Vector{Vec2{Int}}()
    nodeID = Vector{Int}()
    
    for n in 1:nbWires
        helix!(positions, connectivity, nodeID, length, radius*2, nbTotalCells, nbWires, rWireSection, braidingPattern, wireGap, n, true)

    end 
    
    for n in 1:nbWires
        helix!(positions, connectivity, nodeID, length, radius*2, nbTotalCells, nbWires, rWireSection, braidingPattern, wireGap, n, false)
    end


    return positions, connectivity

end

function get_rings(positions)
    
    zvec = Vector{Float64}()
    for i in 1:size(positions,1)
        push!(zvec, positions[i,3])
    end 
    zvec .= round.(zvec, digits = 4)
    unique!(zvec)

    rings = Vector{Vector{Int}}()
    nrings = length(zvec)
    for i in 1:nrings
        ring = Vector{Int}()
        for k in 1:size(positions,1)
            p = positions[k,:]
            if isapprox(round(p[3], digits = 4), zvec[i])
                push!(ring, k)
            end 
        end 
        push!(rings, ring)
    end
    
    return rings 
    
end 

function get_nodespairs_stent(positions)
    
    constraints = Vector{Vec2{Int}}()
    toll = 0.06
    rings = get_rings(positions)

    for r in 1:length(rings)
        
        for i in 1:length(rings[r])
            
            index_i = rings[r][i]
            pos_i = positions[index_i,:]
            
            for j in 1:length(rings[r])
                
                index_j = rings[r][j]
                pos_j = positions[index_j,:]
                
                dist = norm(pos_i-pos_j)
                
                if index_j != index_i && dist < toll && !([index_j, index_i] in constraints)
                    push!(constraints, [index_i, index_j])
                end 
                
            end
        end 
    end 
    
    return constraints
    
end 

# -------------------------------------------------------------------------------------------
# Centerline
# -------------------------------------------------------------------------------------------

function get_centerline_stent(positions_stent, origin=zeros(Vec3))

    rings = get_rings(positions_stent)

    pos_cl_0 = Vector{Vec3{Float64}}()

    nBlocks = length(rings)
    nnodesStent = length(positions_stent)
    indBlock = nnodesStent + 1

    indBlock = nnodesStent + 1
    for ring in 1:nBlocks
        push!(pos_cl_0, [origin[1], origin[2], positions_stent[rings[ring][1]][3]])
        indBlock = indBlock + 1
    end 

    return pos_cl_0
end 

function set_origin!(pos_stent, disp)
    
    pos_stent .-= (disp,)
 
end

function rotate(positions, origin, orient_vessel, orient_stent)
    
    orient_stent = orient_stent/norm(orient_stent)
    orient_vessel = orient_vessel/norm(orient_vessel)
   
    axrot = get_rotation_axis(orient_stent, orient_vessel)
    θ = get_rotation_angle(orient_stent, orient_vessel)
    Mrot = get_rotation_matrix(θ, axrot)

    positions_rot = Vector{Point{3, Float32}}()
    for p in positions
        push!(positions_rot, rotate_around_centre(p, origin, Mrot))
    end 
    
    return positions_rot

end

function get_segments(positions, refp=Vec3(0,0,1))

    npointsCl = length(positions)

    seg = zeros(Vec3, npointsCl)
    segnorm = zeros(Vec3, npointsCl)

    seg[1] = refp

    for i in 2:npointsCl
        seg[i] = positions[i]-positions[i-1]
    end 

    for i in 1:npointsCl
        segnorm[i] = seg[i]/norm(seg[i])
    end 

    return seg, segnorm
    
end 

function get_rotation_axis(a,b)

    if !isapprox(norm(a), 1)
        error("Input vector should be normalised")
    end
    
    if !isapprox(norm(b), 1)
        error("Input vector should be normalised")
    end


    if norm(cross(a,b)) != 0
        return cross(a,b)/norm(cross(a,b))
    else 
        return  [0,0,1]
    end 

end 

function get_rotation_angle(a,b)

    if !isapprox(norm(a), 1)
        error("Input vector should be normalised")
    end
    
    if !isapprox(norm(b), 1)
        error("Input vector should be normalised")
    end
    
    return acos(dot(a,b))
end 

function get_rotation_matrix(θ, axrot)

    c = cos(θ)
    s = sin(θ)
    t = 1-c

    ux = axrot[1]
    uy = axrot[2]
    uz = axrot[3]

    uxx = ux*ux
    uyy = uy*uy
    uzz = uz*uz
    uxy = ux*uy
    uxz = ux*uz
    uyz = uy*uz

    return Mat33(c+uxx*t, uxy*t+uz*s, uxz*t-uy*s, uxy*t-uz*s, c+uyy*t, uyz*t+ux*s, uxz*t+uy*s, uyz*t-ux*s, c+uzz*t)

end 

function rotate_around_centre(p, c, M)

    p = p - c
    p = M * p
    p = p + c

    return p

end 

function write_txt_cl(pos_new, pos_last, step, output_dir)

    if step == 1
        if output_dir != ""
            dir = pwd()
            cd(output_dir)
            foreach(rm, filter(endswith(".vtk"), readdir()))
            foreach(rm, filter(endswith(".txt"), readdir()))
            cd(dir)
        end
    end 

    u = pos_new - pos_last

    open(output_dir * "/u$step.txt", "w") do io
        writedlm(io, u)
    end 

end

# -------------------------------------------------------------------------------------------
# Guides
# -------------------------------------------------------------------------------------------

function get_internal_nodes(positions_stent, rStent, rWireSection, wireGap)
    
    toll = (rWireSection + wireGap/2)/2
    
    Rint = rStent - rWireSection - wireGap/2
    Rext = rStent + rWireSection + wireGap/2
    
    intnodes = []
    extnodes = []
    
    for i in 1:length(positions_stent)
        
        r = sqrt(positions_stent[i][1]^2 + positions_stent[i][2]^2)
        
        if isapprox(r, Rint; atol = toll) 
            push!(intnodes, i)
        elseif isapprox(r, Rext; atol = toll) 
            push!(extnodes, i)
        end 
        
    end 
    
    return intnodes, extnodes
    
end

function get_guide_connectivity_int(rings, nodes_centerline, int_nodes)
    
    if length(rings) != length(nodes_centerline)
        error("length(rings) != length(nodes_centerline).jl")
    end 

    nrings = length(rings)
    constraints = Vector{Vec2{Int}}()
    
    for i in 1:nrings
        thisRing = rings[i]
        nnodesRing = length(thisRing)
        for j in 1:nnodesRing
            if thisRing[j] in int_nodes
                push!(constraints, Vec2(thisRing[j], nodes_centerline[i]))
            end 
        end 
    end     
    
    return constraints
    
end

function build_int_guides_stent_origin(positions, origin, type="Surpass")

    thisWireGap = 0.1

    # total number of stent nodes
    nnodes_stent = length(positions)
    
    # get rings
    rings = get_rings(positions)

    #create centerline 
    positions_guides = Vector{Vec3{Float64}}()
    nguides = length(rings)

    for ring in 1:nguides
        push!(positions_guides, [origin[1], origin[2], positions[rings[ring][1]][3]])
    end
    nnodes_centerline = length(positions_guides)
    nodes_centerline = nnodes_stent+1:nnodes_stent + nnodes_centerline

    # get internal node
    if type == "Surpass"    
        @unpack nbWires, rStent, rWireSection, wireGap, lengthStent, nbTotalCells, braidingPattern = BraidedStent()
        positions_stent, connectivity_stent = compute_bs_geom_given_nbTotalCells(nbWires, rStent, rWireSection, thisWireGap, lengthStent, nbTotalCells, braidingPattern)
        intnodes, extnodes = get_internal_nodes(positions_stent, rStent, rWireSection, thisWireGap*1.1)
    else type == "Wallstent" 
        @unpack nbWires, rStent, rCrimpedStent, rWireSection, lStent, phi, braidingPattern = Wallstent()
        positions_stent, connectivity_stent = compute_bs_geom(nbWires, rStent, rCrimpedStent, rWireSection, thisWireGap, lStent, phi, braidingPattern)
        intnodes, extnodes = get_internal_nodes(positions_stent, rStent, rWireSection, thisWireGap*1.1)
    end

    if length(intnodes) != length(extnodes)
        println("length(intnodes) : ")
        println(length(intnodes))
        println("length(extnodes) : ")
        println(length(extnodes))
        error("Error in the computation of the internal nodes. ")
    end 

    connectivity_guides = get_guide_connectivity_int(rings, nodes_centerline, intnodes)

    return positions_guides, connectivity_guides

end

# -------------------------------------------------------------------------------------------
# Centerline
# -------------------------------------------------------------------------------------------

function read_vtk_centerline(filename)
    
    fid = open(filename, "r")
    
    #skip first four lines of the vtk file 
    readline(fid)
    readline(fid)
    readline(fid)
    readline(fid)
    
    curr_line = split(readline(fid))
    nnodes = parse(Int, curr_line[2])
    
    i = 0
    pts = zeros(nnodes*3)
    
    while i<nnodes*3
        
        curr_line = split(readline(fid))
        n = length(curr_line)
        
        ini = i + 1
        iend = i + n
        index = 1
        
        for j in ini:iend
            pts[j] = parse(Float64, curr_line[index]) 
            index +=1
            
        end 
        
        i = iend
        
    end

    pts_vec = Vector{Vec3}()

    for i in 1:3:nnodes*3

        x = pts[i]
        y = pts[i+1]
        z = pts[i+2]
        push!(pts_vec, Vec3(x,y,z))
    end 

    return unique!(pts_vec)
    
end 

function read_vtk_stent(filename)
    
    fid = open(filename, "r")
    
    #skip first four lines of the vtk file 
    readline(fid)
    readline(fid)
    readline(fid)
    readline(fid)
    
    curr_line = split(readline(fid))
    nnodes = parse(Int, curr_line[2])
    
    i = 0
    pts = zeros(nnodes*3)
    
    while i<nnodes*3
        
        curr_line = split(readline(fid))
        n = length(curr_line)
        
        ini = i + 1
        iend = i + n
        index = 1
        
        for j in ini:iend
            pts[j] = parse(Float64, curr_line[index]) 
            index +=1
            
        end 
        
        i = iend
        
    end

    pts_vec = Vector{Vec3}()

    for i in 1:3:nnodes*3

        x = pts[i]
        y = pts[i+1]
        z = pts[i+2]
        push!(pts_vec, Vec3(x,y,z))
    end 

    return unique!(pts_vec)
    
end 

# -------------------------------------------------------------------------------------------
# Initial conditions
# -------------------------------------------------------------------------------------------

function read_ics_vec(mat)

    vec =  Vector{Vec3{Float64}}()

    for i in 1:size(mat,1)
        push!(vec, (mat[i,:]))
    end 

    return vec

end 

function read_ics_mat(mat)

    vec =  Vector{Mat33{Float64}}()

    for i in 1:size(mat,1)
        push!(vec, (mat[i,:]))
    end 

    return vec

end 

# -------------------------------------------------------------------------------------------
# Write files
# -------------------------------------------------------------------------------------------

function write_txt_centerline(pos_new, pos_last, step, output_dir)

    # if step == 1
    #     if output_dir != ""
    #         dir = pwd()
    #         cd(output_dir)
    #         foreach(rm, filter(endswith(".vtk"), readdir()))
    #         foreach(rm, filter(endswith(".txt"), readdir()))
    #         cd(dir)
    #     end
    # end 

    u = pos_new - pos_last

    open(output_dir * "/u$step.txt", "w") do io
        writedlm(io, u)
    end 

end

function write_vtk_configuration(filename, positions, connectivity)
    
    fid = open(filename, "w")
    
    write(fid,"# vtk DataFile Version 3.0")
    write(fid,"\nvtk output")
    write(fid,"\nASCII")
    write(fid,"\nDATASET POLYDATA")
    
    numberPoints = length(positions)
    write(fid,"\nPOINTS $numberPoints float")  
    for n in positions
        write(fid,"\n")
        write(fid, string.(n[1]))
        write(fid," ")
        write(fid, string.(n[2]))
        write(fid," ")
        write(fid, string.(n[3])) 
    end 
    
    numberLines = length(connectivity)
    numberElementsPerLine = 2
    numberElements= numberLines*(numberElementsPerLine + 1)
    write(fid,"\nLINES $numberLines $numberElements\n")
    
    for i in 1:numberLines
        
        write(fid, string.(numberElementsPerLine))
        write(fid,"\n")
        
        node = connectivity[i][1] -1
        write(fid, string.(node)) 
        write(fid,"\n")
        
        node = connectivity[i][2] -1
        write(fid, string.(node)) 
        write(fid,"\n")
        
    end 
    
    close(fid)
    
end 

function write_vtk_configuration(filename, positions, connectivity, value)
    
    fid = open(filename, "w")
    
    write(fid,"# vtk DataFile Version 3.0")
    write(fid,"\nvtk output")
    write(fid,"\nASCII")
    write(fid,"\nDATASET POLYDATA")
    
    numberPoints = length(positions)
    write(fid,"\nPOINTS $numberPoints float")  
    for n in positions
        write(fid,"\n")
        write(fid, string.(n[1]))
        write(fid," ")
        write(fid, string.(n[2]))
        write(fid," ")
        write(fid, string.(n[3])) 
    end 
    
    numberLines = length(connectivity)
    numberElementsPerLine = 2
    numberElements= numberLines*(numberElementsPerLine + 1)
    write(fid,"\nLINES $numberLines $numberElements\n")
    
    for i in 1:numberLines
        
        write(fid, string.(numberElementsPerLine))
        write(fid,"\n")
        
        node = connectivity[i][1] -1
        write(fid, string.(node)) 
        write(fid,"\n")
        
        node = connectivity[i][2] -1
        write(fid, string.(node)) 
        write(fid,"\n")
        
    end 

    write(fid, "\n POINT_DATA $numberPoints")
    write(fid, "\n SCALARS err float 1")
    write(fid, "\n LOOKUP_TABLE default");
    for i in 1:size(positions, 1)
        write(fid,"\n")
        write(fid, string.(value)) 
    end
    
    close(fid)
    
end 

function write_txt_solution(nodes, beams, nnodes, nbeams, output_dir, CLEAN_FOLDER=true)
    
    if CLEAN_FOLDER
        if output_dir != ""
            dir = pwd()
            cd(output_dir)
            foreach(rm, filter(endswith(".vtk"), readdir()))
            foreach(rm, filter(endswith(".txt"), readdir()))
            foreach(rm, filter(endswith(".vtu"), readdir()))
            foreach(rm, filter(endswith(".txt"), readdir()))
            cd(dir)
        end 
    end

    for i in 1:nnodes
    
        open(output_dir * "/u.txt", "a") do io
            writedlm(io, [nodes.u[i][1] nodes.u[i][2] nodes.u[i][3]])
        end 

        open(output_dir * "/R.txt", "a") do io
            writedlm(io, [nodes.R[i][1,1] nodes.R[i][2,1] nodes.R[i][3,1] nodes.R[i][1,2] nodes.R[i][2,2] nodes.R[i][3,2] nodes.R[i][1,3] nodes.R[i][2,3] nodes.R[i][3,3]] )
        end 

    end 
    
    for i in 1:nbeams

        open(output_dir * "/Re0.txt", "a") do io
            writedlm(io, [beams.Rₑ⁰[i][1,1] beams.Rₑ⁰[i][2,1] beams.Rₑ⁰[i][3,1] beams.Rₑ⁰[i][1,2] beams.Rₑ⁰[i][2,2] beams.Rₑ⁰[i][3,2] beams.Rₑ⁰[i][1,3] beams.Rₑ⁰[i][2,3] beams.Rₑ⁰[i][3,3]])
        end 

    end 
    
end 

function write_vtk_general_structured_mesh(filename, s, dx, x, y, z) 
    
    nx = size(s,1)
    ny = size(s,2)
    nz = size(s,3)
    xo = x[1]
    yo = y[1]
    zo = z[1]
    ncells = nx*ny*nz
    
    # rearrange the sdf in a 1-D vector
    s_vec = zeros(size(s,1)*size(s,2)*size(s,3))
    count = 1
    @inbounds for k in 1:nz
        @inbounds for j in 1:ny
            @inbounds for i in 1:nx
                s_vec[count] = s[i,j,k]
                count = count + 1
            end
        end 
    end 
    
    
    # write vtk
    fid = open(filename, "a")
    
    write(fid,"# vtk DataFile Version 4.2")
    write(fid,"\n vtk output")
    write(fid,"\nASCII")
    write(fid,"\nDATASET STRUCTURED_POINTS")
    write(fid,"\nDIMENSIONS $nx $ny $nz")  
    write(fid,"\nSPACING $dx $dx $dx")  
    write(fid,"\nORIGIN $xo $yo $zo")  
    write(fid,"\nPOINT_DATA  $ncells")  
    write(fid,"\nSCALARS u float 1")  
    write(fid,"\nLOOKUP_TABLE default")  
    write(fid,"\n")  
    
    for n in 1:ncells
        write(fid, string(s_vec[n]))
        write(fid," ")
    end 
    
    close(fid)
    
end 

function get_init_pos_deploy_middle(filename_cl, input, initial_positions_stent, output_dir_crimping)

    positions_stent = initial_positions_stent + read_ics_vec(readdlm(output_dir_crimping * "u.txt"))
    stent = get_centerline_stent(positions_stent)
    η_stent = vcat(0, cumsum(norm.(diff(stent))))
    length_stent = η_stent[end]

    cl = read_vtk_centerline(filename_cl)
    η = vcat(0, cumsum(norm.(diff(cl))))
    length_cl = η[end]

    up_bound = length_cl/2
    down_bound = length_cl/2-length_stent
    deploy_pos = down_bound + input*(up_bound-down_bound)

    return deploy_pos

end 

function save_obj(filename, vertices, trisconn)

    fid = open(filename, "w")

    for p in vertices
        write(fid, "v ", string(p[1]), " ", string(p[2]), " ", string(p[3]))
        write(fid,"\n")
    end


    for f in trisconn
        write(fid, "f ", join(convert.(Int, f), " "))
        write(fid,"\n")
    end

    close(fid)

end