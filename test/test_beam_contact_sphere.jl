function test_sphere()
    
    # -------------------------------------------------------------------------------------------
    # Building the nodes
    # -------------------------------------------------------------------------------------------
    
    # positions
    L = 10/100; dx = L/10
    
    X₀ =  Vec3{T}[]
    
    for x in 0:dx:L
        push!(X₀, Vec3(0, 0, x))
    end
    
    # total number of nodes
    nnodes = size(X₀,1)
    
    # initial conditions
    u⁰ = zeros(size(X₀,1)*3)
    u̇⁰ = zeros(size(X₀,1)*3)
    ü⁰ = zeros(size(X₀,1)*3)
    w⁰ = zeros(size(X₀,1)*3)
    ẇ⁰ = zeros(size(X₀,1)*3)
    ẅ⁰ = zeros(size(X₀,1)*3)
    
    # plane for cylindrical coordinates
    plane = fill("xy", length(X₀))
    
    # nodes StructArray
    nodes = constructor_nodes(X₀, u⁰, u̇⁰, ü⁰, w⁰, ẇ⁰, ẅ⁰, nothing, T)
    
    # -------------------------------------------------------------------------------------------
    # Building the beams
    # -------------------------------------------------------------------------------------------
    
    # total number of beams
    nbeams = nnodes #!!!
    
    # conn
    conn = Vec2{Int}[]
    aux1 = 1:nbeams-1
    aux2 = 2:nbeams
    
    for i in 1:nbeams-1
        push!(conn, (aux1[i], aux2[i])) 
    end
    
    # interpolation points per beam
    nbInterpolationPoints = 30
    
    # geometric dimension of the cylinder
    d = 0.6/1000
    r = d/2
    E = 5*1e7
    nu = 0.33
    ρ = 7850
    G = E/(2*(1+nu))
    A = pi*r^2
    I₂₂ = 0.25*pi*r^4
    I₃₃ = 0.25*pi*r^4
    Iₒ = I₂₂+I₃₃
    Iᵣᵣ = Iₒ
    J = Iₒ
    Jᵨ = Mat33(ρ*Iₒ, 0, 0, 0, ρ*I₂₂, 0, 0, 0, ρ*I₃₃)
    Aᵨ = ρ*A
    
    geom = Geometry{T}(A, I₂₂, I₃₃, Iₒ, Iᵣᵣ, J)
    mat = Material{T}(E, G, Aᵨ, Jᵨ)
    
    # beams vector
    beams = constructor_beams(nodes, conn, mat, geom, nbInterpolationPoints, nothing)
    
    #-----------------------------------------------------------------------------------
    # Simulation parameters
    # -------------------------------------------------------------------------------------------
    
    # integration parameters
    α = -0.05
    β = 0.25*(1-α)^2
    γ = 0.5*(1-2*α)
    damping = 0
    
    # time step and total time
    Δt = 0.01
    Δt_plot = 0.1
    tᵉⁿᵈ =  2.5
    
    # tolerance and maximum number of iterations
    tol_res = 1e-5
    tol_ΔD = 1e-5
    max_it = 10
    
    # Gauss points
    nG = 3
    ωG = Vec3(5/9, 8/9, 5/9)
    zG = Vec3(-sqrt(3/5), 0, sqrt(3/5)) 
    
    # penalty parameters
    kₙ = 50
    μ = 0.3
    εᵗ = 0.1
    
    comp = SimulationParameters(α, β, γ, damping,  Δt, Δt_plot, tᵉⁿᵈ, tol_res, tol_ΔD, max_it, nG, ωG, zG, kₙ, μ, εᵗ, T)
    
    # -------------------------------------------------------------------------------------------
    # External forces
    # -------------------------------------------------------------------------------------------
    
    # external force and applied dof
    flag_crimping = false
    F(t) = 0
    dofs_load = T[]
    
    ext_forces = ExternalForces(flag_crimping, F, dofs_load, T)
    
    # -------------------------------------------------------------------------------------------
    # Boundary conditions
    # -------------------------------------------------------------------------------------------
    
    # multifreedom constraints
    cons = nothing
    
    # number of dof (6 per node)
    ndofs = nnodes*6
    
    # Dirichlet boundary conditions: fixed positions
    fixed_dofs = 1:6
    free_dofs = setdiff(1:ndofs, fixed_dofs) 
    
    # Dirichlet boundary conditions: moving positions
    flag_cylindrical = false
    disp_dofs = [1]
    dispA = 10/100
    tA = 2.5
    k = dispA/tA
    Fdisp(t) = (k*t).*(t.<=tA) .+ (-k*t .+ 2*dispA).*((t.>tA) .& (t.<=2*tA)) .+ (k*(t.-2*tA)).*((t.>2*tA) .& (t.<=3*tA)) .+ (-k*t .+ 4*dispA).*((t.>3*tA) .& (t.<=4*tA))
    flag_disp_vector = false
    disp_vals = T[]
    
    # boundary conditions strucutre 
    bcs = BoundaryConditions(fixed_dofs, free_dofs, flag_cylindrical, flag_disp_vector, Fdisp, disp_vals, disp_dofs, T)
    
    # -------------------------------------------------------------------------------------------
    # SDF
    # -------------------------------------------------------------------------------------------
    
    # sphere center
    x0 = 0.05
    y0 = 0.005
    z0 = 0.08
    
    # sphere radius
    R = 0.04
    
    sdf = Sphere_SDF{T}(r, R, x0, y0, z0)
    
    # -------------------------------------------------------------------------------------------
    # Configuration
    # -------------------------------------------------------------------------------------------
    
    # configuration: mesh, external forces and boundary conditions
    conf = Configuration(mat, geom, nnodes, ndofs, ext_forces, bcs, T)
    
    # -------------------------------------------------------------------------------------------
    # Solve
    # -------------------------------------------------------------------------------------------
    
    params = ParamsTest()
    solver!(nodes, beams, conf, comp, sdf, cons, params, T)         
    
    # -------------------------------------------------------------------------------------------
    # Test 
    # -------------------------------------------------------------------------------------------
    
    rtol = 1e-2
    
    pos_matlab = Vec3(0.107853165032904, 0.0132224390067658, 0.0360871249073232)
    @test isapprox(nodes.X₀[5] + nodes.u[5], pos_matlab; rtol)
    
    pos_matlab = Vec3(0.129725721716302, 0.0548254754528272, 0.0730574700480992)
    @test isapprox(nodes.X₀[end] + nodes.u[end], pos_matlab; rtol)
    
end 

test_sphere()