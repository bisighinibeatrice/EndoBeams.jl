function test_angle()

    #-------------------------------------------------------------------------------------------
    # Building the nodes
    # -------------------------------------------------------------------------------------------

    # positions
    dx = 2.5; L = 10

    X₀ =  Vec3{T}[]
    for x in 0:dx:L
        push!(X₀, Vec3(0, x, 0))
    end
    for x in dx:dx:L
        push!(X₀, Vec3(-x,  L, 0))
    end

    # total number of nodes
    nnodes = size(X₀,1)

    # initial conditions
    u⁰ = zeros(nnodes*3)
    u̇⁰ = zeros(nnodes*3)
    ü⁰ = zeros(nnodes*3)
    w⁰ = zeros(nnodes*3)
    ẇ⁰ = zeros(nnodes*3)
    ẅ⁰ = zeros(nnodes*3)

    # plane for cylindrical coordinates
    plane = fill("xy", length(X₀))

    # nodes StructArray
    nodes = constructor_nodes(X₀, u⁰, u̇⁰, ü⁰, w⁰, ẇ⁰, ẅ⁰, plane, nothing, T)

    # -------------------------------------------------------------------------------------------
    # Building the beams
    # -------------------------------------------------------------------------------------------

    # total number of beams
    nbeams = nnodes-1

    # conn
    conn = Vec2{Int}[]
    aux1 =  1:nnodes-1
    aux2 =  2:nnodes

    for i in 1:nbeams
        push!(conn, (aux1[i], aux2[i])) 
    end

    # interpolation points per beam
    nbInterpolationPoints = 30

    # geometric and material properties
    E = 1e6
    G = 1e6
    Jᵨ = Mat33(20, 0, 0, 0, 10, 0, 0, 0, 10)
    Aᵨ = 1
    A = 1
    I₂₂ = 1e-3
    I₃₃ = 1e-3
    Iₒ = 0
    Iᵣᵣ = 0
    J = 1e-3

    geom = Geometry{T}(A, I₂₂, I₃₃, Iₒ, Iᵣᵣ, J)
    mat = Material{T}(E, G, Aᵨ, Jᵨ)

    # beams vector
    allbeams = constructor_beams(nodes, conn, mat, geom, nbInterpolationPoints, nothing)

    #-----------------------------------------------------------------------------------
    # Simulation parameters
    # -------------------------------------------------------------------------------------------

    # integration parameters
    α = -0.05
    β = 0.25*(1-α)^2
    γ = 0.5*(1-2*α)
    damping = 0

    # time step and total time
    dt = 0.25
    dt_plot = 0.25
    tend = 150

    # tolerance and maximum number of iterations
    tol_res = 1e-5
    tol_ddk = 1e-5
    max_it = 10

    # Gauss points
    nG = 3
    ωG = Vec3(5/9, 8/9, 5/9)
    zG = Vec3(-sqrt(3/5), 0, sqrt(3/5)) 

    # penalty parameters
    eps_C = 5000
    μ = 0
    εₜ = 0.1

    comp = constructor_simulation_parameters(α, β, γ, damping,  dt, dt_plot, tend, tol_res, tol_ddk, max_it, nG, ωG, zG, eps_C, μ, εₜ, T)

    # -------------------------------------------------------------------------------------------
    # External forces
    # -------------------------------------------------------------------------------------------

    # external force and applied dof
    flag_crimping = false
    Fext(t) = 1*(t.<=1).*(50*t) .+1*((t.>1) .& (t.<=2)).*(-50*t.+100).+(t.>2).*0
    dofs_load = Int[]
    push!(dofs_load, 6*size((0:dx:L),1)-3)

    ext_forces = constructor_ext_force(flag_crimping, Fext, dofs_load, T)

    # -------------------------------------------------------------------------------------------
    # Boundary conditions
    # -------------------------------------------------------------------------------------------

    # multifreedom constrains
    cons = T[]

    # Dirichlet boundary conditions: fixed positions
    ndofs = nnodes*6
    fixed_dofs = 1:6
    free_dofs = setdiff(1:ndofs, fixed_dofs)

    # Dirichlet boundary conditions: moving positions
    flag_cylindrical = false
    Fdisp(t) = 0
    dofs_disp = Int[]
    flag_disp_vector = false
    udisp = T[]

    # boundary conditions strucutre 
    bcs = constructor_boundary_conditions(fixed_dofs, free_dofs, flag_cylindrical, flag_disp_vector, Fdisp, udisp, dofs_disp, T)

    # -------------------------------------------------------------------------------------------
    # SDF
    # -------------------------------------------------------------------------------------------

    #there is no contact
    sdf = nothing

    # -------------------------------------------------------------------------------------------
    # Configuration
    # -------------------------------------------------------------------------------------------

    # configuration: mesh, external forces and boundary conditions
    conf = constructor_configuration(mat, geom, nnodes, ndofs, ext_forces, bcs, T)

    # -------------------------------------------------------------------------------------------
    # Solve
    # -------------------------------------------------------------------------------------------

    params = ParamsTest()
    solver!(nodes, allbeams, conf, comp, sdf, cons, params, T)       

    # -------------------------------------------------------------------------------------------
    # Test
    # -------------------------------------------------------------------------------------------
  
    pos_matlab = Vec3(-10.6867358894953, 6.61144877793072, -6.2584504990589)
    @test isapprox(nodes.X₀[end] + nodes.u[end], pos_matlab; atol=1e-6)

end

test_angle()