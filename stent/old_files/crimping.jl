function crimping(rStent, rCrimpedStent, initial_positions_stent, connectivity_stent, output_dir_crimping="output3D/outputCrimping3D")

    # -------------------------------------------------------------------------------------------
    # Time stepping parameters
    # -------------------------------------------------------------------------------------------

    # initial time step and total time
    ini_Δt = 0.1
    max_Δt = 1.
    Δt_plot =  0.1
    tᵉⁿᵈ = 1

    params = Params(;ini_Δt, Δt_plot, max_Δt, tᵉⁿᵈ, output_dir = output_dir_crimping, stop_on_energy_threshold=true, energy_threshold=1e-12, tol_res = 1e-3, tol_ΔD = 1e-3, record_timings=false, verbose=true)

    # -----------------------------------------------------------------------------------------
    # Building the nodes
    # -------------------------------------------------------------------------------------------

    # total number of nodes
    nnodes = size(initial_positions_stent, 1)

    # initial conditions
    u⁰ = zeros(Vec3, nnodes)
    u̇⁰ = zeros(Vec3, nnodes)
    ü⁰ = zeros(Vec3, nnodes)
    w⁰ = zeros(Vec3, nnodes)
    ẇ⁰ = zeros(Vec3, nnodes)
    ẅ⁰ = zeros(Vec3, nnodes)

    plane = "xy"

    # nodes StructArray
    nodes = build_nodes(initial_positions_stent, u⁰, u̇⁰, ü⁰, w⁰, ẇ⁰, ẅ⁰, plane, nothing)

    # -------------------------------------------------------------------------------------------
    # Building the beams
    # -------------------------------------------------------------------------------------------

    # geometric and material properties
    E = 225*1e3
    ν = 0.33
    mass_scaling = 1e3
    ρ = 9.13*1e-9 * mass_scaling
    radius = 0.014
    damping = 1e4

    # beams vector
    beams = build_beams(nodes, connectivity_stent, E, ν, ρ, radius, damping, nothing)
    nbeams = length(beams)

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
    nodespairs = get_nodespairs_stent(initial_positions_stent)
    constraints = build_constraints(nodespairs, kᶜᵒⁿ, ηᶜᵒⁿ)

    # Dirichlet boundary conditions: blocked positions
    aux1 = 1:6:ndofs-5
    aux2 = 2:6:ndofs-4
    fixed_dofs = sort(unique([aux1; aux2]))
    free_dofs = setdiff(1:ndofs, fixed_dofs)

    # Dirichlet dof (x6)
    use_cylindrical_coords = true
    disp_dofs = 1:6:ndofs-5
    disp_vals = zeros(ndofs)
    dispA = -(rStent-rCrimpedStent)
    tA = tᵉⁿᵈ
    k = dispA/tA
    disp(t, node_idx) =  (k.*t).*(t.<tA) .+ (k.*tA).*(t.>=tA)

    # boundary conditions strucutre
    bcs = BoundaryConditions(fixed_dofs, free_dofs, disp, disp_vals, disp_dofs, use_cylindrical_coords)

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
    write_txt_solution(nodes, beams, nnodes, nbeams, output_dir_crimping)
    
end 

