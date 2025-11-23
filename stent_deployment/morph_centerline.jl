"""
Morph a crimped stent centerline toward its deployed configuration.

This function takes a crimped stent geometry (from a previous crimping simulation) 
and iteratively morphs its centerline so that it aligns with a target (deployed) 
centerline loaded from a VTK file. The morphing is performed via incremental 
rotations of centerline segments, controlled by a specified number of iterations 
and deployment offset.

# Arguments
- `free_positions::Matrix{Float64}`: Initial nodal positions of the stent.
- `n_iterations::Int`: Number of morphing iterations (controls smoothness of deformation).
- `deployment_offset_mm::Float64`: Offset along the target centerline (mm).
- `target_centerline_filename::String`: Path to deployed centerline VTK file.
- `output_dir_crimping::String`: Directory containing crimping simulation results.
- `output_dir_morphing_centerline::String`: Directory to store morphing results.

# Returns
- `deployment_origin_point::Vec3`: Reference position along the stent centerline indicating where deployment starts (first node of the crimped centerline).
"""
function morph_stent_centerline_to_deployment!(
    free_positions,
    n_iterations,
    deployment_offset_mm,
    target_centerline_filename,
    output_dir_crimping,
    output_dir_morphing_centerline
)

    #-------------------------------
    # LOAD TARGET CENTERLINE
    #-------------------------------
    target_centerline = read_vtk_centerline_points(target_centerline_filename)

    # Parameterize target centerline by arc length
    η_target = vcat(0, cumsum(norm.(diff(target_centerline))))
    spline_target = Vec3(Dierckx.Spline1D(η_target, getindex.(target_centerline, i)) for i in 1:3)

    #-------------------------------
    # RECONSTRUCT CRIMPED STENT CONFIGURATION
    #-------------------------------
    crimped_positions = matrix_to_vec3_array(free_positions) +
                        matrix_to_vec3_array(readdlm(output_dir_crimping * "u.txt"))

    crimped_centerline = compute_stent_centerline(crimped_positions)

    n_centerline_nodes = length(crimped_centerline)
    connectivity_centerline = [Vec2(i, i + 1) for i in 1:n_centerline_nodes - 1]

    #-------------------------------
    # BUILD TARGET MORPHED CENTERLINE
    #-------------------------------
    η_initial = vcat(0, cumsum(norm.(diff(crimped_centerline))))
    target_morphed_centerline = zeros(Vec3, n_centerline_nodes)

    for i in 1:n_centerline_nodes
        target_morphed_centerline[i] = Dierckx.evaluate.(spline_target, deployment_offset_mm + η_initial[i])
    end

    #-------------------------------
    # ALIGN INITIAL CONFIGURATION TO TARGET
    #-------------------------------
    shift_positions!(crimped_positions, crimped_centerline[1] - target_morphed_centerline[1])
    crimped_centerline = compute_stent_centerline(crimped_positions, target_morphed_centerline[1])

    #-------------------------------
    # SEGMENT DIRECTIONS & ROTATIONS
    #-------------------------------
    target_segments, target_directions = compute_segments_and_directions(target_morphed_centerline)

    rotation_axes   = [compute_rotation_axis(target_directions[i], target_directions[i - 1]) for i in 2:n_centerline_nodes]
    rotation_angles = [compute_rotation_angle(target_directions[i], target_directions[i - 1]) for i in 2:n_centerline_nodes]

    max_rotations   = -rotation_angles
    delta_rotations = max_rotations ./ n_iterations
    current_rotations = zeros(Float64, length(rotation_angles))

    current_nodes = copy(crimped_centerline)
    morphed_nodes = copy(crimped_centerline)

    # Save initial morphing state
    write_vtk_polyline(output_dir_morphing_centerline * "/centerline_step_0.vtk", morphed_nodes, connectivity_centerline)
    writedlm(output_dir_morphing_centerline * "/u_step_0.txt", morphed_nodes .- crimped_centerline)

    #-------------------------------
    # ITERATIVE MORPHING LOOP
    #-------------------------------
    for iter = 1:n_iterations
        current_nodes .= crimped_centerline
        current_rotations .+= delta_rotations

        # Apply incremental rotations segment by segment
        for seg = 1:length(current_rotations)
            rotation_matrix = rodrigues_rotation_matrix(current_rotations[seg], rotation_axes[seg])

            # Keep upstream nodes fixed
            for j = 1:seg
                morphed_nodes[j] = current_nodes[j]
            end

            # Rotate downstream nodes
            for j = seg + 1:length(current_nodes)
                morphed_nodes[j] = rotate_point_around_center(current_nodes[j], current_nodes[seg], rotation_matrix)
            end

            current_nodes .= morphed_nodes
        end

        # Save current morphing step
        write_vtk_polyline(output_dir_morphing_centerline * "/centerline_step_$iter.vtk", morphed_nodes, connectivity_centerline)
        writedlm(output_dir_morphing_centerline * "/u_step_$iter.txt", morphed_nodes .- crimped_centerline)
    end

    # Return reference point of stent deployment (first node of crimped centerline)
    deployment_origin_point = crimped_centerline[1]

    return deployment_origin_point
end
