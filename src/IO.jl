
#----------------------------------
# WRITE VTK FILES
#----------------------------------

# Write nodes info as vtk file
function write_VTK_nodes(step, allnodes, allbeams, dirOutput) 
    
    filename=string(dirOutput * "/nodes$step.vtk");
     
    fid = open(filename, "a")

    write(fid,"# vtk DataFile Version 3.0")
    write(fid,"\n vtk output")
    write(fid,"\nASCII")
    write(fid,"\nDATASET POLYDATA")
    
    numberPoints = length(allnodes)
    write(fid,"\nPOINTS $numberPoints float")  
    for n in allnodes
        write(fid,"\n")
        write(fid, string.(n.pos[1] + n.u[1])) 
        write(fid," ")
        write(fid, string.(n.pos[2] + n.u[2])) 
        write(fid," ")
        write(fid, string.(n.pos[3] + n.u[3])) 
    end 
    
    numberLines = length(allbeams)
    numberElementsPerLine = 2
    numberElements= numberLines*(numberElementsPerLine + 1)
    write(fid,"\nLINES $numberLines $numberElements\n")
    
    for i in 1:numberLines

        write(fid, string.(numberElementsPerLine))
        write(fid,"\n")

        node = allbeams.node1[i] -1
        write(fid, string.(node)) 
        write(fid,"\n")

        node = allbeams.node2[i] -1
        write(fid, string.(node)) 
        write(fid,"\n")

    end 

    write(fid, "\n POINT_DATA $numberPoints")
    write(fid, "\n SCALARS u float 3")
    write(fid, "\n LOOKUP_TABLE default");
    for n in allnodes
        write(fid,"\n")
        write(fid, string.(n.u[1])) 
        write(fid," ")
        write(fid, string.(n.u[2])) 
        write(fid," ")
        write(fid, string.(n.u[3])) 
    end

    close(fid)
    
end 

# Write interpolated beam info as vtk file
function write_VTK_beams(step, allnodes, allbeams, positions, connectivity, dirOutput) 

    get_centreline!(positions, connectivity, allnodes, allbeams)
    
    filename=string(dirOutput * "/beams$step.vtk");
     
    fid = open(filename, "a")

    write(fid,"# vtk DataFile Version 3.0")
    write(fid,"\nvtk output")
    write(fid,"\nASCII")
    write(fid,"\nDATASET POLYDATA")
    
    numberPoints = size(positions, 1)
    write(fid,"\nPOINTS $numberPoints float")  
    for i in 1: numberPoints
        write(fid,"\n")
        write(fid, string.(positions[i,1])) 
        write(fid," ")
        write(fid, string.(positions[i,2])) 
        write(fid," ")
        write(fid, string.(positions[i,3])) 
    end 
    
    numberLines = size(connectivity, 2)
    numberElementsPerLine= size(connectivity, 1)
    numberElements= numberLines*(numberElementsPerLine+1)
    write(fid,"\nLINES $numberLines $numberElements\n")
    
    for i in 1:numberLines

        write(fid, string.(numberElementsPerLine))
        write(fid,"\n")

        for j in 1:numberElementsPerLine
            node::Int64 = connectivity[j,i] -1
            write(fid, string.(node)) 
            write(fid,"\n")
        end

    end 

    close(fid)
    
end 

function get_connectivity_centerline(pos_cl)

    conn = Vector{Vec2{Int}}()
    aux1 = 1:size(pos_cl,1)-1
    aux2 = 2:size(pos_cl,1)
    for i in 1:size(pos_cl,1)-1
        push!(conn, (aux1[i], aux2[i])) 
    end

    return conn

end 

# Write Gaussian points info as vtk file
function write_VTK_GP(step, sol_GP, dirOutput) 
    
    filename = string(dirOutput * "/GP$step.vtk");
    fid = open(filename, "a")

    write(fid,"# vtk DataFile Version 3.0")
    write(fid,"\nvtk output")
    write(fid,"\nASCII")
    write(fid,"\nDATASET POLYDATA")
    
    numberPoints = length(sol_GP.xGP)
    write(fid,"\nPOINTS $numberPoints float")  
    for i in 1: numberPoints
        write(fid,"\n")
        write(fid, string.(sol_GP.xGP[i][1])) 
        write(fid," ")
        write(fid, string.(sol_GP.xGP[i][2])) 
        write(fid," ")
        write(fid, string.(sol_GP.xGP[i][3])) 
    end 
    
    conn = get_connectivity_centerline(sol_GP.xGP)

    numberLines = size(conn, 1)
    numberElementsPerLine = 2
    numberElements= numberLines*(numberElementsPerLine+1)
    write(fid,"\nLINES $numberLines $numberElements\n")
    
    for i in 1:numberLines

        write(fid, string.(numberElementsPerLine))
        write(fid,"\n")
        node1::Int64 = conn[i][1] - 1
        write(fid, string.(node1)) 
        write(fid,"\n")
        node2::Int64 = conn[i][2] - 1
        write(fid, string.(node2)) 
        write(fid,"\n")
    end

    write(fid,"\nPOINT_DATA $numberPoints")
    write(fid,"\nSCALARS fN float 3")
    write(fid,"\nLOOKUP_TABLE default")
    for i in 1: numberPoints
        write(fid,"\n")
        write(fid, string.(sol_GP.fGP_N[i][1])) 
        write(fid,"\n")
        write(fid, string.(sol_GP.fGP_N[i][2]))
        write(fid,"\n") 
        write(fid, string.(sol_GP.fGP_N[i][3])) 
    end 

    write(fid,"\nSCALARS fT float 3")
    write(fid,"\nLOOKUP_TABLE default")
    write(fid,"\n")
    for i in 1: numberPoints
        write(fid,"\n")
        write(fid, string.(sol_GP.fGP_T[i][1])) 
        write(fid,"\n")
        write(fid, string.(sol_GP.fGP_T[i][2]))
        write(fid,"\n") 
        write(fid, string.(sol_GP.fGP_T[i][3])) 
    end

    write(fid,"\nSCALARS g float 1")
    write(fid,"\nLOOKUP_TABLE default")
    for i in 1: numberPoints
        write(fid,"\n")
        write(fid, string.(sol_GP.gGP[i][1])) 
    end

    write(fid,"\nSCALARS status float 1")
    write(fid,"\nLOOKUP_TABLE default")
    for i in 1: numberPoints
        write(fid,"\n")
        write(fid, string.(sol_GP.status[i][1])) 
    end

    close(fid)
    
end 

# Write nodes as vtk file (no allnodes, no allbeams)
function write_VTK_configuration(filename, X, C)
    
    fid = open(filename, "a")
    
    write(fid,"# vtk DataFile Version 3.0")
    write(fid,"\nvtk output")
    write(fid,"\nASCII")
    write(fid,"\nDATASET POLYDATA")
    
    numberPoints = length(X)
    write(fid,"\nPOINTS $numberPoints float")  
    for n in X
        write(fid,"\n")
        write(fid, string.(n[1]))
        write(fid," ")
        write(fid, string.(n[2]))
        write(fid," ")
        write(fid, string.(n[3])) 
    end 
    
    numberLines = length(C)
    numberElementsPerLine = 2
    numberElements= numberLines*(numberElementsPerLine + 1)
    write(fid,"\nLINES $numberLines $numberElements\n")
    
    for i in 1:numberLines
        
        write(fid, string.(numberElementsPerLine))
        write(fid,"\n")
        
        node = C[i][1] -1
        write(fid, string.(node)) 
        write(fid,"\n")
        
        node = C[i][2] -1
        write(fid, string.(node)) 
        write(fid,"\n")
        
    end 
    
    close(fid)
    
end 

#----------------------------------
# READ TXT  FILES
#----------------------------------

"""
pos = read_TXT_file_pos(filename)
Read a .txt file as a Vector{Vec3{T}}() for position and displacement. 
"""
function read_TXT_file_pos(filename)

    Xpc = Vector{Vec3{T}}()
    xpc = readdlm(filename)

    for i in 1:size(xpc,1)
        push!(Xpc, Vec3{T}(xpc[i,1], xpc[i,2], xpc[i,3]))
    end

    return Xpc

end 

"""
conn = read_TXT_file_pos(filename)
Read a .txt file as a Vector{Vec2{T}}() for connectivity. 
"""
function read_TXT_file_conn(filename)

    Xpc = Vector{Vec2{Int}}()
    xpc = readdlm(filename)

    for i in 1:size(xpc,1)
        if xpc[i,1] < xpc[i,2]
            push!(Xpc, Vec2{Int}(xpc[i,1], xpc[i,2]))
        else 
            push!(Xpc, Vec2{Int}(xpc[i,2], xpc[i,1]))
        end 
    end

    return Xpc

end 

"""
Ics_vec = read_TXT_file_ICs_array(filename)
Read a .txt file as a Vector{T}() for displacement, velocity and acceleration initial conditions. 
"""
function read_TXT_file_ICs_array(filename)

    X = readdlm(filename)
    nnodes3 = size(X,1)*size(X,2)
    x = zeros(nnodes3)

    j = 1
    for i in 1:3:nnodes3
        x[i] = X[j,1]
        x[i+1] = X[j,2]
        x[i+2] = X[j,3]
        j = j+1
    end

    return x

end

"""
Ics_mat = read_TXT_file_ICs_matrix(filename)
Read a .txt file as a Vector{T}() for rotation matrix (nodes and edges) initial conditions. 
"""
function read_TXT_file_ICs_matrix(filename)

    X = readdlm(filename)
    nnodes = size(X,1)
    x = Vector{Mat33{T}}()

    for i in 1:nnodes
       push!(x, Mat33{T}(X[i,:]))
    end

    return x

end

#----------------------------------
# READ TXT  FILES
#----------------------------------

# Read a discrete SDF from a VTK file
function read_VTK_sdf(filename)
    
    fid = open(filename, "r")
    
    #skip first four lines of the vtk file 
    readline(fid)
    readline(fid)
    readline(fid)
    readline(fid)
    
    curr_line = split(readline(fid))
    
    nnodex = parse(Int, curr_line[2]) 
    nnodey = parse(Int, curr_line[3]) 
    nnodez = parse(Int, curr_line[4]) 
    
    curr_line = split(readline(fid))
    
    dx = parse(Float64, curr_line[2]) 
    dy = parse(Float64, curr_line[3]) 
    dz = parse(Float64, curr_line[4]) 
    
    curr_line = split(readline(fid))
    
    x0 = parse(Float64, curr_line[2]) 
    y0 = parse(Float64, curr_line[3]) 
    z0 = parse(Float64, curr_line[4]) 
    xend = x0 + (nnodex-1)*dx
    yend = y0 + (nnodey-1)*dy
    zend = z0 + (nnodez-1)*dz
    dom = [x0, xend, y0, yend, z0, zend]
    
    curr_line = split(readline(fid))
    nnode = parse(Int, curr_line[2])
    
    #skip other two lines
    readline(fid)
    readline(fid)
    
    i = 0
    sdf = zeros(nnode)
    
    while i<nnode
        
        curr_line = split(readline(fid))
        n = length(curr_line)
        
        ini  = i + 1
        iend = i + n
        index = 1
        
        for j in ini:iend
            sdf[j] = parse(Float64, curr_line[index]) 
            index +=1
            
        end 
        
        i = iend
        
    end
    
    return nnodex, nnodey, nnodez, dx, dy, dz, dom, sdf
    
end 
