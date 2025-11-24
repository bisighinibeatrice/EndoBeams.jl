# -------------------------------------------------------------------------------------------
# Convert variables
# -------------------------------------------------------------------------------------------

"""
Convert a vector of `Vec3{Float64}` to a matrix of shape `n × 3`.

# Arguments
- `arr::Vector{Vec3{Float64}}`: Vector of 3D points.

# Returns
- `Matrix{Float64}`: Matrix with 3 columns (x, y, z coordinates), where each row 
  corresponds to a `Vec3` in the input vector.
"""
function vec3_array_to_matrix(arr)
    n = length(arr)
    mat = Matrix{Float64}(undef, n, 3)
    for i in 1:n
        v = arr[i]
        mat[i, 1] = v.x
        mat[i, 2] = v.y
        mat[i, 3] = v.z
    end
    return mat
end

"""
Convert a matrix of shape `n × 3` to a vector of `Vec3{Float64}`.

# Arguments
- `mat::Matrix{Float64}`: Matrix with 3 columns (x, y, z coordinates).

# Returns
- `Vector{Vec3{Float64}}`: Each row of the matrix is converted to a `Vec3`.
"""
function matrix_to_vec3_array(mat)
    @assert size(mat, 2) == 3 "Input matrix must have exactly 3 columns."

    return [Vec3(mat[i, 1], mat[i, 2], mat[i, 3]) for i in 1:size(mat, 1)]
end

"""
Convert a matrix of shape `n × 2` to a vector of `Vec2{Int}`.

# Arguments
- `mat::Matrix{Int}`: Matrix with 2 columns.

# Returns
- `Vector{Vec2{Int}}`: Each row of the matrix is converted to a `Vec2`.
"""
function matrix_to_vec2_array(mat)
    @assert size(mat, 2) == 2 "Input matrix must have exactly 2 columns."

    return [Vec2(mat[i, 1], mat[i, 2]) for i in 1:size(mat, 1)]
end

"""
Convert a matrix of shape `n × 9` to a vector of `Mat33{Float64}` matrices.

# Arguments
- `mat::Matrix{Float64}`: Matrix with 9 columns representing 3×3 matrices row-wise.

# Returns
- `Vector{Mat33{Float64}}`: Each row of the matrix converted to a `Mat33`.
"""
function matrix_to_mat33_array(mat::Matrix{Float64})
    @assert size(mat, 2) == 9 "Input matrix must have exactly 9 columns."

    vec = Vector{Mat33{Float64}}(undef, size(mat, 1))

    for i in 1:size(mat, 1)
        # # reshape row vector into 3x3 matrix in column-major order
        # # note: Julia is column-major, so transpose if needed
        # # here, we assume row-major in input, so reshape and transpose
        # row = mat[i, :]
        # mat33 = reshape(row, 3, 3)'  # transpose to get proper orientation
        vec[i] = Mat33(mat[i, :])
    end

    return vec
end

# -------------------------------------------------------------------------------------------
# Read
# -------------------------------------------------------------------------------------------

"""
Read a `.vtk` centerline file and return the node coordinates as unique Vec3 points.

# Arguments
- `filename::String`: Path to the `.vtk` file.

# Returns
- `points::Vector{Vec3}`: Unique list of 3D points along the centerline.
"""
function read_vtk_centerline_points(filename::String)
    fid = open(filename, "r")

    # Skip the VTK header (first 4 lines)
    for _ in 1:4
        readline(fid)
    end

    # Parse number of nodes from the next line (format: "POINTS <n> float")
    point_info = split(readline(fid))
    n_nodes = parse(Int, point_info[2])

    coords = Float64[]
    total_coords = 3 * n_nodes

    # Read all coordinates from the file
    while length(coords) < total_coords
        line_values = split(readline(fid))
        append!(coords, parse.(Float64, line_values))
    end

    # Reshape into Vec3 format
    points = Vector{Vec3{Float64}}()
    for i in 1:3:total_coords
        push!(points, Vec3(coords[i], coords[i + 1], coords[i + 2]))
    end

    close(fid)
    return unique!(points)
end

# -------------------------------------------------------------------------------------------
# Write
# -------------------------------------------------------------------------------------------

"""
Write node and beam solution data (`u`, `R`, `Re0`) to text files.

Arguments:
- `nodes`: Object containing displacement `u` and rotation matrices `R` (per node).
- `beams`: Object containing reference element rotations `Rₑ⁰` (per beam).
- `num_nodes`: Number of nodes.
- `num_beams`: Number of beams.
- `output_dir`: Directory where files will be written.
- `CLEAN_FOLDER`: If true, clears existing `.vtk`, `.vtu`, `.txt` files in the folder.
"""
function export_stent_solution_to_txt(nodes, beams, num_nodes, num_beams, output_dir::String, flag_clean_folder::Bool=true)
    
    # Clean existing output files if requested
    if flag_clean_folder && output_dir != ""
        savedir = pwd()
        cd(output_dir)
        for ext in [".vtk", ".vtu", ".txt"]
            for file in filter(f -> endswith(f, ext), readdir())
                rm(file, force=true)
            end
        end
        cd(savedir)
    end

    # Ensure directory path ends without double slash
    if !endswith(output_dir, "/")
        output_dir *= "/"
    end

    # Write node displacements and rotations
    for i in 1:num_nodes
        # Write displacement vector u (3 values)
        open(output_dir * "u.txt", "a") do io
            writedlm(io, [nodes.u[i]'])  # Transpose to row
        end

        # Write rotation matrix R (flattened 3x3)
        R = nodes.R[i]
        open(output_dir * "R.txt", "a") do io
            writedlm(io, [vec(R)'])  # Flatten and transpose to row
        end
    end

    # Write beam reference rotations
    for i in 1:num_beams
        Re0 = beams.Rₑ⁰[i]
        open(output_dir * "Re0.txt", "a") do io
            writedlm(io, [vec(Re0)'])  # Flatten 3x3 matrix to row
        end
    end
end

"""
Write a `.vtk` file in ASCII POLYDATA format from node positions and connectivity.

# Arguments
- `filename::String`: Output file path.
- `positions::Vector{Vec3}`: List of 3D point positions.
- `connectivity::Vector{Vec2{Int}}`: Each pair defines a line segment (1-based indexing).
"""
function write_vtk_polyline(filename::String, positions, connectivity)

    open(filename, "w") do fid
        # VTK header
        write(fid, "# vtk DataFile Version 3.0\n")
        write(fid, "vtk output\n")
        write(fid, "ASCII\n")
        write(fid, "DATASET POLYDATA\n")

        # Write point coordinates
        n_points = length(positions)
        write(fid, "POINTS $n_points float\n")
        for p in positions
            write(fid, "$(p[1]) $(p[2]) $(p[3])\n")
        end

        # Write line connectivity
        n_lines = length(connectivity)
        n_cells = n_lines * 3  # 3 = 2 point indices + 1 count per line
        write(fid, "LINES $n_lines $n_cells\n")

        for line in connectivity
            # VTK uses 0-based indexing
            a, b = line .- 1
            write(fid, "2 $a $b\n")
        end
    end
end
