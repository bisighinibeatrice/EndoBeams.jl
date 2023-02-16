function positioning(initial_positions_stent, connectivity_stent, nb_iterations, origin, output_dir_crimping, output_dir_positioning_cl, output_dir_positioning)
    
    # -------------------------------------------------------------------------------------------
    # Time stepping parameters
    # -------------------------------------------------------------------------------------------
    
    # initial time step and total time
    ini_Δt = 1.
    max_Δt = 1.
    Δt_plot =  0.0
    tᵉⁿᵈ = nb_iterations-1
    
    params = Params(;ini_Δt, Δt_plot, max_Δt, tᵉⁿᵈ, output_dir = output_dir_positioning, stop_on_energy_threshold=true, energy_threshold=1e-12, tol_res = 1e-3, tol_ΔD = 1e-3, record_timings=false, verbose=false)
    
    # -----------------------------------------------------------------------------------------
    # Stent 
    # -------------------------------------------------------------------------------------------
    
    crimped_positions_stent = initial_positions_stent .+  read_ics_vec(readdlm(output_dir_crimping * "u.txt"))
    positions_cl =  get_centerline_stent(crimped_positions_stent)
    set_origin!(crimped_positions_stent, positions_cl[1]-origin)
    positions_cl =  get_centerline_stent(crimped_positions_stent, origin)

    # -----------------------------------------------------------------------------------------
    # Guides 
    # -------------------------------------------------------------------------------------------
    
    positions_guides, connectivity_guides = build_int_guides_stent_origin(crimped_positions_stent, positions_cl[1])

    # -----------------------------------------------------------------------------------------
    # Guides = stent
    # -------------------------------------------------------------------------------------------
    
    initial_positions = [crimped_positions_stent; positions_guides]
    connectivity = [connectivity_stent; connectivity_guides]

    nnodes_stent = length(crimped_positions_stent)
    nbeams_stent = length(connectivity_stent)
    nnodes = length(initial_positions)  
    nbeams = length(connectivity)  
    iguides = nnodes_stent+1:nnodes

    # -----------------------------------------------------------------------------------------
    # Building the nodes 
    # -------------------------------------------------------------------------------------------
    
    # total number of nodes
    nnodes = size(initial_positions, 1)
    
    # initial conditions
    u⁰ = zeros(Vec3, nnodes)
    u̇⁰ = zeros(Vec3, nnodes)
    ü⁰ = zeros(Vec3, nnodes)
    w⁰ = zeros(Vec3, nnodes)
    ẇ⁰ = zeros(Vec3, nnodes)
    ẅ⁰ = zeros(Vec3, nnodes)
    
    plane = "xy"

    R⁰ =  zeros(Mat33, nnodes)
    R⁰[1:nnodes_stent] .= read_ics_mat(readdlm(output_dir_crimping * "R.txt"))
    R⁰[nnodes_stent+1:nnodes] .= (ID3,)
    
    # nodes StructArray
    nodes = build_nodes(initial_positions, u⁰, u̇⁰, ü⁰, w⁰, ẇ⁰, ẅ⁰, plane, R⁰)
    
    # -------------------------------------------------------------------------------------------
    # Building the beams
    # -------------------------------------------------------------------------------------------
    
    # geometric and material properties
    E = 225*1e3
    ν = 0.33
    ρ = 9.13*1e-3
    radius = 0.065
    damping = 1E2

    # E = 225*1e3
    # ν = 0.33
    # mass_scaling = 1E6
    # ρ = 9.13*1e-9 * mass_scaling
    # radius = 0.014
    # damping = 1E3*7.5

    Re₀ =  zeros(Mat33, nbeams)
    Re₀[1:nbeams_stent] .= read_ics_mat(readdlm(output_dir_crimping * "Re0.txt"))
    for (i,c) in enumerate(connectivity_guides)
        node1 = nodes[c[1]]
        node2 = nodes[c[2]]
        Re₀[nbeams_stent+i] = get_Rₑ⁰(node1.X₀, node2.X₀)
    end 

    # beams vector
    beams = build_beams(nodes, connectivity, E, ν, ρ, radius, damping, Re₀)
    
    # contact parameters
    kₙ = 4/3 * 5/(1-0.5^2)*sqrt(radius) # Approximate Hertz contact with 5 MPa wall stiffness
    μ = 0.01
    εᵗ = 0.001 #regularized parameter for friction contact
    ηₙ = 0.1
    kₜ = kₙ
    ηₜ = ηₙ
    u̇ₛ = 0.01

    contact = ContactParameters(kₙ, μ, εᵗ, ηₙ, kₜ, ηₜ, u̇ₛ)
 
    # -------------------------------------------------------------------------------------------
    # External forces
    # -------------------------------------------------------------------------------------------
    
    # external force and applied dof
    loaded_dofs = []
    force(t, node_idx) = 0
    
    ext_forces = ExternalForces(force, loaded_dofs)
    
    # -------------------------------------------------------------------------------------------
    # Boundary conditions
    # -------------------------------------------------------------------------------------------
    
    # number of dof (6 per node)
    ndofs = nnodes*6
    
    # penalty constraints
    kᶜᵒⁿ = 1e3
    ηᶜᵒⁿ = 1
    nodespairs = get_nodespairs_stent(initial_positions)
    constraints = build_constraints(nodespairs, kᶜᵒⁿ, ηᶜᵒⁿ)
    
    # Dirichlet boundary conditions: blocked positions
    fixed_dofs = Vector{Int}() 
    for i in iguides
        push!(fixed_dofs, 6*(i-1)+1)
        push!(fixed_dofs, 6*(i-1)+2)
        push!(fixed_dofs, 6*(i-1)+3)
    end
    free_dofs = setdiff(1:ndofs, fixed_dofs)
    
    # Dirichlet dof (x6)
    flag_cylindrical = false
    flag_load_from_file = true
    dir_folder_load = output_dir_positioning_cl
    disp_dofs = fixed_dofs
    disp_vals = zeros(ndofs)
    disp(t, node_idx) =  0

    # boundary conditions strucutre
    bcs = BoundaryConditions(fixed_dofs, free_dofs, disp, disp_vals, disp_dofs, flag_cylindrical, flag_load_from_file, dir_folder_load)
    
    # -------------------------------------------------------------------------------------------
    # SDF
    # -------------------------------------------------------------------------------------------
    
    sdf = nothing
    
    # -------------------------------------------------------------------------------------------
    # Final configuration
    # -------------------------------------------------------------------------------------------
    
    # configuration: mesh, external forces and boundary conditions
    conf = Configuration(nodes, beams, constraints, ext_forces, bcs, contact, sdf)
    
    # -------------------------------------------------------------------------------------------
    # Start simulation
    # -------------------------------------------------------------------------------------------
    
    solver!(conf, params)
    write_txt_solution(nodes, beams, nnodes_stent, nbeams_stent, output_dir_positioning, false)

end 
