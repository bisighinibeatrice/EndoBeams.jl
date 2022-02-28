

# Write interpolated beam info as vtk file
function write_VTK(write_counter, step, t, nodes, beams, energy, conf, sdf, comp, vtkdata) 

    fill!(vtkdata.normal_contact_force, zero(eltype(vtkdata.normal_contact_force)))
    fill!(vtkdata.tangential_contact_force, zero(eltype(vtkdata.tangential_contact_force)))
    fill!(vtkdata.incontact, 0)

    recompute_at_gausspts!(vtkdata, nodes, beams, conf.mat, sdf, comp)

    output_dir = vtkdata.output_dir

    vtk_grid("$output_dir/data_$write_counter", vtkdata.interpolated_points, vtkdata.interpolated_lines; compress = false) do vtk

        vtk["Stress", VTKPointData()] = vtkdata.stress
        vtk["Strain", VTKPointData()] = vtkdata.strain

        if !isnothing(sdf)
            vtk["Contact gap", VTKPointData()] = vtkdata.contact_distance
            vtk["Normal contact force", VTKPointData()] = vtkdata.normal_contact_force
            vtk["Tangential contact force", VTKPointData()] = vtkdata.tangential_contact_force
            vtk["In contact?", VTKPointData()] = vtkdata.incontact
        end

        vtk["Displacement", VTKPointData()] = vtkdata.displacement
        vtk["Velocity", VTKPointData()] = vtkdata.velocity

        vtk["Contact Energy", VTKFieldData()] = energy.contact_energy
        vtk["Kinetic Energy", VTKFieldData()] = energy.kinetic_energy
        vtk["Strain Energy", VTKFieldData()] = energy.strain_energy

        # vtk["time", VTKFieldData()] = t
        # vtk["step", VTKFieldData()] = step

        vtkdata.VTKcollection[t] = vtk 

    end

     
    
end 


function recompute_at_gausspts!(vtkdata, nodes, beams, mat, sdf, comp)

    @batch for bi in eachindex(beams)

        n1 = beams.node1[bi]
        n2 = beams.node2[bi]

        X₁, X₂ = nodes.X₀[n1], nodes.X₀[n2]
        u₁, u₂ = nodes.u[n1], nodes.u[n2]
        u̇₁, u̇₂ = nodes.u̇[n1], nodes.u̇[n2]
        ẇ₁, ẇ₂ = nodes.ẇ[n1], nodes.ẇ[n2]
        R₁, R₂ = nodes.R[n1], nodes.R[n2]
        l₀ = beams.l₀[bi]
        Rₑ⁰ = beams.Rₑ⁰[bi]


        x₁ =  X₁ + u₁
        x₂ =  X₂ + u₂
        
        lₙ = norm(x₂ - x₁)

        Rₑ, _, _, _, Gᵀ¹, Gᵀ², Gᵀ³, Gᵀ⁴, _ = local_Rₑ_and_aux(x₁, x₂, R₁, R₂, Rₑ⁰[:,2], lₙ)

        R̅₁ = Rₑ' * R₁ * Rₑ⁰
        R̅₂ = Rₑ' * R₂ * Rₑ⁰

        Θ̅₁ = toangle(R̅₁)
        Θ̅₂ = toangle(R̅₂)

        
        @inbounds for (i, zᴳ) in enumerate(range(-1, 1, length=vtkdata.intermediate_points))

            k = vtkdata.intermediate_points * (bi-1) + i 

            ξ = l₀*(zᴳ+1)/2

            # Shape functions
            N₁ = 1-ξ/l₀
            N₂ = 1-N₁
            N₃ = ξ*(1-ξ/l₀)^2
            N₄ = -(1-ξ/l₀)*((ξ^2)/l₀)

            uᵗ = @SVector [0, N₃*Θ̅₁[3] + N₄*Θ̅₂[3], -N₃*Θ̅₁[2] + -N₄*Θ̅₂[2]]
            xᴳ = N₁*x₁ + N₂*x₂ + Rₑ*uᵗ

            vtkdata.interpolated_points[k] = xᴳ

            vtkdata.displacement[k] = N₁*u₁ + N₂*u₂
            vtkdata.velocity[k] = N₁*u̇₁ + N₂*u̇₂

            vtkdata.strain[k] = 1-lₙ/l₀
            vtkdata.stress[k] = mat.E*(1-lₙ/l₀)

            if !isnothing(sdf)
                pₙ, _, _, gₙ, ∂gₙ∂x, _ =  contact_gap(xᴳ, comp.εᶜ, sdf)
                vtkdata.contact_distance[k] = gₙ
                if pₙ > 0 
                    U̇₁ = Rₑ' * u̇₁
                    U̇₂ = Rₑ' * u̇₂
                    Ẇ₁ = Rₑ' * ẇ₁
                    Ẇ₂ = Rₑ' * ẇ₂
                    Suᵗ = skew(uᵗ)
                    N₇ = N₃+N₄
                    N₇lₙ = N₇/lₙ
                    P₁P¹ = @SMatrix [0 0 0; 0 N₇lₙ 0;0 0 N₇lₙ]
                    P₁P² = @SMatrix [0 0 0; 0 0 N₃;0 -N₃ 0]
                    P₁P³ = -P₁P¹
                    P₁P⁴ = @SMatrix [0 0 0; 0 0 N₄;0 -N₄ 0]
                    H₁¹ = N₁*ID3 + P₁P¹ - Suᵗ*Gᵀ¹
                    H₁² =          P₁P² - Suᵗ*Gᵀ²
                    H₁³ = N₂*ID3 + P₁P³ - Suᵗ*Gᵀ³
                    H₁⁴ =          P₁P⁴ - Suᵗ*Gᵀ⁴
                    h₁ = H₁¹ * U̇₁ + H₁² * Ẇ₁ + H₁³ * U̇₂ + H₁⁴ * Ẇ₂
                    u̇₀ = Rₑ * h₁
                    ġₙ = dot(u̇₀, ∂gₙ∂x)*∂gₙ∂x
                    ġₜ = u̇₀ - ġₙ
                    ġₜ² = dot(ġₜ, ġₜ)
                    γᵈᵃᵐᵖ = comp.γᵈᵃᵐᵖ
                    μʳᵉᵍ = comp.μ/sqrt(ġₜ²+comp.εᵗ)
                    vtkdata.normal_contact_force[k] = pₙ*∂gₙ∂x - γᵈᵃᵐᵖ * pₙ * ġₙ
                    vtkdata.tangential_contact_force[k] = -μʳᵉᵍ * pₙ * ġₜ
                    vtkdata.incontact[k] = 1
                end
            end

            k += 1

        end

    end

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
    dom = (x0, xend, y0, yend, z0, zend)
    
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


