# Structure to hold data required for VTK output visualization
struct VTKDataBeams
    VTKcollection::WriteVTK.CollectionFile # Collection of VTK files
    output_dir::String                    # Directory for VTK output files
    intermediate_points::Int              # Number of intermediate points in beam interpolation
    interpolated_points::Vector{Vec3{Float64}} # Interpolated points along the beams
    interpolated_lines::Vector{MeshCell{WriteVTK.PolyData.Lines, UnitRange{Int}}} # Interpolated lines along the beams
    stress::Vector{Float64}               # Stress values at each interpolated point
    strain::Vector{Float64}               # Strain values at each interpolated point
    displacement::Vector{Vec3{Float64}}   # Displacement at each interpolated point
    velocity::Vector{Vec3{Float64}}       # Velocity at each interpolated point
    contact_distance::Vector{Float64}     # Contact distance at each interpolated point
    normal_contact_force::Vector{Vec3{Float64}} # Normal contact force at each interpolated point
    tangential_contact_force::Vector{Vec3{Float64}} # Tangential contact force at each interpolated point
    tangential_velocity::Vector{Float64}  # Tangential velocity at each interpolated point
    incontact::Vector{Int}                # Contact status indicator
    status::Vector{Int}                   # Contact status (active or inactive)

    nodes::Vector{Vec3{Float64}} 
    lines_nodes::Vector{MeshCell{WriteVTK.PolyData.Lines, Vector{Int}}} 
    displacement_nodes::Vector{Vec3{Float64}} 
    incontact_nodes::Vector{Int}                
    contactforce_nodes::Vector{Vec3{Float64}} 
    contactdistance_nodes::Vector{Float64}        
end

# Constructor for VTKData structure, initializing interpolation and contact data
function VTKDataBeams(conf, output_dir, sdf = Nothing, intermediate_points = 30)

    num_beams = length(conf.beams)
    num_nodes = length(conf.nodes)

    # Initialize vectors for interpolation and contact data
    interpolated_points = zeros(Vec3{Float64}, num_beams * intermediate_points)
    interpolated_lines = [MeshCell(PolyData.Lines(), (i - 1) * intermediate_points + 1 : i * intermediate_points) for i in 1:num_beams]
    stress = zeros(length(interpolated_points))
    strain = zeros(length(interpolated_points))
    displacement = zeros(Vec3{Float64}, length(interpolated_points))
    velocity = zeros(Vec3{Float64}, length(interpolated_points))

    # Conditional initialization of contact-related vectors if SignedDistanceField is provided
    if !isnothing(sdf)
        contact_distance = zeros(length(interpolated_points))
        normal_contact_force = zeros(Vec3{Float64}, length(interpolated_points))
        tangential_contact_force = zeros(Vec3{Float64}, length(interpolated_points))
        tangential_velocity = zeros(length(interpolated_points))
        incontact = zeros(Int, length(interpolated_points))
        status = zeros(Int, length(interpolated_points))
    else
        contact_distance = Float64[]
        normal_contact_force = Vec3{Float64}[]
        tangential_contact_force = Vec3{Float64}[]
        tangential_velocity = Float64[]
        incontact = Int[]
        status = Int[]
    end

    # Initialize vectors for interpolation and contact data
    nodes = zeros(Vec3{Float64}, num_nodes)
    lines_nodes = [MeshCell(PolyData.Lines(), [beam.node1, beam.node2]) for beam in conf.beams]
    displacement_nodes = zeros(Vec3{Float64}, num_nodes)
    incontact_nodes = zeros(Int, num_nodes)              
    contactforce_nodes = zeros(Vec3{Float64}, num_nodes)
    contactdistance_nodes = zeros(Float64, num_nodes)
    
    # Prepare output directory and initialize Paraview collection
    clean_folders(output_dir)
    collection = paraview_collection("$output_dir/simulation")

    return VTKDataBeams(collection, output_dir, intermediate_points, interpolated_points, interpolated_lines, stress, strain, displacement, velocity, contact_distance, normal_contact_force, tangential_contact_force, tangential_velocity, incontact, status, nodes, lines_nodes, displacement_nodes, incontact_nodes, contactforce_nodes, contactdistance_nodes)
end

# Writes VTK files for beams simulation results at a given timestep.
function write_VTK_beams(write_counter, step, t, conf, vtkdata)
    # Unpack configuration
    @unpack nodes, beams = conf

    # Reset contact-related data
    recompute_at_gausspts_beams!(vtkdata, nodes, beams)

    # Output directory for VTK files
    output_dir = vtkdata.output_dir

    # Write the VTK data for this timestep
    vtk_grid("$output_dir/beams_$write_counter", vtkdata.interpolated_points, vtkdata.interpolated_lines; compress = false) do vtk
      
        # Write point-based data (stress, strain, etc.)
        vtk["Stress", VTKPointData()] = vtkdata.stress
        vtk["Strain", VTKPointData()] = vtkdata.strain

        # Write other results (displacement, velocity, energy, etc.)
        vtk["Displacement", VTKPointData()] = vtkdata.displacement
        vtk["Velocity", VTKPointData()] = vtkdata.velocity
        # vtk["Kinetic EnergySimulationState", VTKFieldData()] = energyⁿ⁺¹.kinetic_energy
        # vtk["Strain EnergySimulationState", VTKFieldData()] = energyⁿ⁺¹.strain_energy

        # Time and step information
        vtk["time", VTKFieldData()] = t
        vtk["step", VTKFieldData()] = step

        # Store the VTK data for this time step
        vtkdata.VTKcollection[t] = vtk
    end

    # Reset contact-related data
    recompute_at_nodes_beams!(vtkdata, conf.nodes)

    # Write the VTK data for this timestep
    vtk_grid("$output_dir/beams_nodes_$write_counter", vtkdata.nodes, vtkdata.lines_nodes; compress = false) do vtk
       
        vtk["Displacement", VTKPointData()] = vtkdata.displacement_nodes
        vtk["In Contact", VTKPointData()] = vtkdata.incontact_nodes
        vtk["Contact Force", VTKPointData()] = vtkdata.contactforce_nodes
        vtk["Contact Distance", VTKPointData()] = vtkdata.contactdistance_nodes

        vtk["time", VTKFieldData()] = t
        vtk["step", VTKFieldData()] = step

        # Store the VTK data for this time step
        vtkdata.VTKcollection[t] = vtk
    end
end

# Function to recompute results at Gauss points for each beam element
function recompute_at_gausspts_beams!(vtkdata, nodes, beams)
    # Loop over each beam element
    @batch for bi in eachindex(beams)
        # Unpack node and beam data
        n1 = beams.node1[bi]
        n2 = beams.node2[bi]
        X₁, X₂ = nodes.X₀[n1], nodes.X₀[n2]
        u₁, u₂ = nodes.u[n1], nodes.u[n2]
        u̇₁, u̇₂ = nodes.u̇[n1], nodes.u̇[n2]
        ẇ₁, ẇ₂ = nodes.ẇ[n1], nodes.ẇ[n2]
        R₁, R₂ = nodes.R[n1], nodes.R[n2]
        l₀ = beams.l₀[bi]
        Rₑ⁰ = beams.Rₑ⁰[bi]

        # Compute beam geometry and updated positions
        x₁ = X₁ + u₁
        x₂ = X₂ + u₂
        lₙ = norm(x₂ - x₁)

        # Compute rotation matrices and other auxiliary variables
        Rₑ, _, _, _, Gᵀ¹, Gᵀ², Gᵀ³, Gᵀ⁴, _ = local_Rₑ_and_aux(x₁, x₂, R₁, R₂, Rₑ⁰[:,2], lₙ)
        R̅₁ = Rₑ' * R₁ * Rₑ⁰
        R̅₂ = Rₑ' * R₂ * Rₑ⁰
        Θ̅₁ = toangle(R̅₁)
        Θ̅₂ = toangle(R̅₂)

        # Loop over Gauss points and compute quantities
        @inbounds for (i, zᴳ) in enumerate(range(-1, 1, length=vtkdata.intermediate_points))
            k = vtkdata.intermediate_points * (bi - 1) + i

            # Shape function interpolation
            ξ = l₀ * (zᴳ + 1) / 2
            N₁ = 1 - ξ / l₀
            N₂ = 1 - N₁
            N₃ = ξ * (1 - ξ / l₀)^2
            N₄ = -(1 - ξ / l₀) * (ξ^2 / l₀)

            # Displacement and velocity interpolation
            uᵗ = @SVector [0, N₃ * Θ̅₁[3] + N₄ * Θ̅₂[3], -N₃ * Θ̅₁[2] - N₄ * Θ̅₂[2]]
            xᴳ = N₁ * x₁ + N₂ * x₂ + Rₑ * uᵗ

            # Store interpolated results
            vtkdata.interpolated_points[k] = xᴳ
            vtkdata.displacement[k] = N₁ * u₁ + N₂ * u₂
            vtkdata.velocity[k] = N₁ * u̇₁ + N₂ * u̇₂
            vtkdata.strain[k] = 1 - lₙ / l₀
            vtkdata.stress[k] = beams[bi].properties.E * (1 - lₙ / l₀)

        end
    end
end

# Recomputes results at Gauss points for each solid element.
function recompute_at_nodes_beams!(vtkdata, nodes)
    @batch for n in eachindex(nodes)
        vtkdata.nodes[n] = nodes.X₀[n] + nodes.u[n]  # Updated position
        vtkdata.displacement_nodes[n] = nodes.u[n]  
        vtkdata.incontact_nodes[n] = nodes.incontact[n]  
        vtkdata.contactforce_nodes[n] = nodes.contactforce[n]  
        vtkdata.contactdistance_nodes[n] = nodes.contactdistance[n]  
    end 
end
