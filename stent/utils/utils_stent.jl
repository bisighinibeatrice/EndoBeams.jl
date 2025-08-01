
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

# -------------------------------------------------------------------------------------------
# Centerline
# -------------------------------------------------------------------------------------------

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

# -------------------------------------------------------------------------------------------
# Centerline
# -------------------------------------------------------------------------------------------

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
# Write files
# -------------------------------------------------------------------------------------------


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

function get_init_pos_deploy_middle(filename_cl, input, initial_positions, output_dir_crimping)

    positions = initial_positions + matrix_to_vec3_array(readdlm(output_dir_crimping * "u.txt"))
    stent = compute_stent_centerline(positions)
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