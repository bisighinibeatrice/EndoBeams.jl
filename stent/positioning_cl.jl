function positioning_cl(initial_positions_stent, connectivity_stent, nb_iterations, deploy_pos, filename_cl, output_dir_crimping, output_dir_positioning_cl)

    dir = pwd()
    cd(output_dir_positioning_cl)
    foreach(rm, filter(endswith(".txt"), readdir()))
    foreach(rm, filter(endswith(".vtk"), readdir())) 
    cd(dir)

    #-------------------------------------------------------------

    # read final centerline
    # cl_T = reverse(read_vtk_centerline(filename_cl))
    cl_T = read_vtk_centerline(filename_cl)

    v_η = vcat(0, cumsum(norm.(diff(cl_T))))
    spline = Vec3(Dierckx.Spline1D(v_η, getindex.(cl_T, i)) for i in 1:3)

    # get initial centerline
    positions = initial_positions_stent + read_ics_vec(readdlm(output_dir_crimping * "u.txt"))
    pos_cl_0 = get_centerline_stent(positions)

    nnodesCl = length(pos_cl_0)
    aux1 = 1:nnodesCl-1
    aux2 = 2:nnodesCl

    connectivity_cl = Vec2{Int}[]
    for i in 1:nnodesCl-1
        push!(connectivity_cl, (aux1[i], aux2[i])) 
    end

    # get final centerline
    cl_η = vcat(0, cumsum(norm.(diff(pos_cl_0))))
    npointsCl = length(cl_η)
    pos_cl_T = zeros(Vec3, npointsCl)
    for i in 1:npointsCl
        pos_cl_T[i] = Dierckx.evaluate.(spline, deploy_pos + cl_η[i])
    end 
    write_vtk_configuration(output_dir_positioning_cl * "/pos_cl_T.vtk", pos_cl_T, connectivity_cl)

    # centre stent in deploy_pos
    set_origin!(positions, pos_cl_0[1]-pos_cl_T[1])   
    write_vtk_configuration(output_dir_positioning_cl * "/crimped_stent.vtk", positions, connectivity_stent)

    pos_cl_0 =  get_centerline_stent(positions, pos_cl_T[1])
    write_vtk_configuration(output_dir_positioning_cl * "/pos_cl_0.vtk", pos_cl_0, connectivity_cl)

    #-------------------------------------------------------------

    seg_T, segnorm_T = get_segments(pos_cl_T)
    #-------------------------------------------------------------

    axrot_T = zeros(Vec3, npointsCl-1)
    for i in 2:npointsCl
       axrot_T[i-1] = get_rotation_axis(segnorm_T[i], segnorm_T[i-1])
    end 

    #-------------------------------------------------------------

    θ_T = zeros(Float64, npointsCl-1)
    for i in 2:npointsCl
        θ_T[i-1] = get_rotation_angle(segnorm_T[i], segnorm_T[i-1])
    end 

    #-------------------------------------------------------------

    θmax = -θ_T
    dθ = θmax/nb_iterations
    θnew = 0*θ_T
    nodes_old = similar(pos_cl_0)
    nodes_i =  zeros(Vec3, length(nodes_old))

    for n = 1:nb_iterations  

        nodes_old .= pos_cl_0
        θnew = θnew + dθ

        for i in 1:length(θnew)
            
            Mrot_i = get_rotation_matrix(θnew[i], axrot_T[i])

            for j in 1:i
                nodes_i[j] = nodes_old[j]
            end 
            
            for k in (i+1):length(nodes_old)
                p = rotate_around_centre(nodes_old[k], nodes_old[i], Mrot_i)
                nodes_i[k] = p
            end 
            
            nodes_old .= nodes_i

        end 

        write_vtk_configuration(output_dir_positioning_cl * "/clMorph_$n.vtk", nodes_i, connectivity_cl)
        write_txt_centerline(nodes_i, pos_cl_0, n, output_dir_positioning_cl)
        
    end

    return pos_cl_0
end