function read_vtk_tetrahedral_mesh(file_path::String)

    # Open the file for reading
    open(file_path, "r") do f
        # Initialize storage for nodes and connectivity
        nodes = []
        tetrahedra_connectivity = []
        triangle_connectivity = []
        
        # Flags to know which section we're reading
        in_points_section = false
        in_cells_section = false
        num_points = 0
        num_cells = 0

        # Read through the lines of the file
        for line in eachline(f)

            # Skip empty lines or comments
            if isempty(line) || startswith(line, "#")
                continue
            end

            # Check if the line contains the keyword "POINTS" indicating the start of the points section
            if occursin("POINTS", line)
                in_points_section = true
                num_points = parse(Int, split(line)[2])  # Get the number of points
                continue
            end

            # Read nodes if we are in the points section
            if in_points_section
                if  length(nodes) < num_points
                    node_coords = parse.(Float64, split(line))
                    push!(nodes, node_coords)
                end
                if length(nodes) == num_points
                    in_points_section = false  # End of POINTS section
                end
                continue
            end

            # Check if the line contains the keyword "CELLS" indicating the start of the cells section
            if occursin("CELLS", line)
                in_cells_section = true
                num_cells = parse(Int, split(line)[2])  # Get the number of cells
                continue
            end

            # Read cell data if we are in the cells section
            if in_cells_section
                # Process the line containing cell data
                cell_data = parse.(Int, split(line))
                num_cell_points = cell_data[1]  # The first number indicates the number of points in the cell
                cell_connectivity = cell_data[2:end] .+ 1  # The remaining numbers are the node indices for the cell

                # Store the connectivity
                if num_cell_points == 4
                    # Tetrahedron connectivity (4 nodes)
                    push!(tetrahedra_connectivity, cell_connectivity)
                elseif num_cell_points == 3
                    # Triangle connectivity (3 nodes)
                    push!(triangle_connectivity, cell_connectivity)
                end

                # If we've read all cells, we can stop processing the cells section
                if length(tetrahedra_connectivity) + length(triangle_connectivity) == num_cells
                    break
                end
            end
        end
        
        nodes  = vcat(nodes'...)  
        tetrahedra_connectivity  = vcat(tetrahedra_connectivity'...)  
        triangle_connectivity  = vcat(triangle_connectivity'...)  

        return nodes, tetrahedra_connectivity, triangle_connectivity
    end
end

function read_vtk_triangle_mesh(file_path::String)
    # Open the file for reading
    open(file_path, "r") do f
        # Initialize storage for nodes and connectivity
        nodes = []
        triangle_connectivity = []

        # Flags to know which section we're reading
        in_points_section = false
        in_cells_section = false
        num_points = 0
        num_cells = 0

        # Read through the lines of the file
        for line in eachline(f)

            # Skip empty lines or comments
            if isempty(line) || startswith(line, "#")
                continue
            end

            # Check if the line contains the keyword "POINTS" indicating the start of the points section
            if occursin("POINTS", line)
                in_points_section = true
                num_points = parse(Int, split(line)[2])  # Get the number of points
                continue
            end

            # Read nodes if we are in the points section
            if in_points_section
                if length(nodes) < num_points
                    node_coords = parse.(Float64, split(line))
                    push!(nodes, node_coords)
                end
                if length(nodes) == num_points
                    in_points_section = false  # End of POINTS section
                end
                continue
            end

            # Check if the line contains the keyword "CELLS" indicating the start of the cells section
            if occursin("CELLS", line)
                in_cells_section = true
                num_cells = parse(Int, split(line)[2])  # Get the number of cells
                continue
            end

            # Read cell data if we are in the cells section
            if in_cells_section
                # Process the line containing cell data
                cell_data = parse.(Int, split(line))
                num_cell_points = cell_data[1]  # The first number indicates the number of points in the cell
                cell_connectivity = cell_data[2:end] .+ 1  # The remaining numbers are the node indices for the cell

                # Store the connectivity for triangles (3 nodes)
                if num_cell_points == 3
                    push!(triangle_connectivity, cell_connectivity)
                end

                # If we've read all cells, we can stop processing the cells section
                if length(triangle_connectivity) == num_cells
                    break
                end
            end
        end
        
        nodes  = vcat(nodes'...)  
        triangle_connectivity  = vcat(triangle_connectivity'...)  

        return  [(nodes[i, :]) for i in 1:size(nodes, 1)],  [(triangle_connectivity[i, :]) for i in 1:size(triangle_connectivity, 1)]
    end
end
