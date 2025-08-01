function morph_stent_centerline_to_deployment!(
    crimped_positions_matrix,
    stent_connectivity,
    n_iterations,
    deployment_offset,
    target_centerline_filename,
    crimping_results_dir,
    output_dir
)

    # Load target (deployed) centerline from file
    target_centerline = read_vtk_centerline_points(target_centerline_filename)

    η_target = vcat(0, cumsum(norm.(diff(target_centerline))))
    spline_target = Vec3(Dierckx.Spline1D(η_target, getindex.(target_centerline, i)) for i in 1:3)

    # Reconstruct crimped stent configuration
    crimped_positions = matrix_to_vec3_array(crimped_positions_matrix) +
                        matrix_to_vec3_array(readdlm(crimping_results_dir * "u.txt"))

    initial_centerline = compute_stent_centerline(crimped_positions)

    n_centerline_nodes = length(initial_centerline)
    connectivity_cl = [Vec2(i, i + 1) for i in 1:n_centerline_nodes - 1]

    # Generate target centerline based on deployment offset
    η_initial = vcat(0, cumsum(norm.(diff(initial_centerline))))
    target_morphed_centerline = zeros(Vec3, n_centerline_nodes)

    for i in 1:n_centerline_nodes
        target_morphed_centerline[i] = Dierckx.evaluate.(spline_target, deployment_offset + η_initial[i])
    end

    write_vtk_polyline(output_dir * "/target_centerline.vtk", target_morphed_centerline, connectivity_cl)

    # Align crimped stent to start at deployed position
    shift_positions!(crimped_positions, initial_centerline[1] - target_morphed_centerline[1])
    write_vtk_polyline(output_dir * "/aligned_crimped_stent.vtk", crimped_positions, matrix_to_vec2_array(stent_connectivity))

    initial_centerline = compute_stent_centerline(crimped_positions, target_morphed_centerline[1])
    write_vtk_polyline(output_dir * "/centerline_step_0.vtk", initial_centerline, connectivity_cl)

    # Compute segment directions for the target centerline
    target_segments, target_directions = compute_segments_and_directions(target_morphed_centerline)

    # Compute rotation axes and angles between segment directions
    rotation_axes = [compute_rotation_axis(target_directions[i], target_directions[i - 1]) for i in 2:n_centerline_nodes]
    rotation_angles = [compute_rotation_angle(target_directions[i], target_directions[i - 1]) for i in 2:n_centerline_nodes]

    max_rotations = -rotation_angles
    delta_rotations = max_rotations ./ n_iterations
    current_rotations = zeros(Float64, length(rotation_angles))

    current_nodes = copy(initial_centerline)
    morphed_nodes = similar(initial_centerline)

    # Iteratively morph the centerline
    for iter = 1:n_iterations
        current_nodes .= initial_centerline
        current_rotations .+= delta_rotations

        for seg = 1:length(current_rotations)
            rotation_matrix = rodrigues_rotation_matrix(current_rotations[seg], rotation_axes[seg])

            for j = 1:seg
                morphed_nodes[j] = current_nodes[j]
            end

            for j = seg + 1:length(current_nodes)
                morphed_nodes[j] = rotate_point_around_center(current_nodes[j], current_nodes[seg], rotation_matrix)
            end

            current_nodes .= morphed_nodes
        end

        write_vtk_polyline(output_dir * "/centerline_step_$iter.vtk", morphed_nodes, connectivity_cl)
        writedlm(output_dir * "/u_step_$iter.txt", morphed_nodes .- initial_centerline)
    end

    return initial_centerline
end
