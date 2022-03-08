function test_ring_crimping()
    
    # -------------------------------------------------------------------------------------------
    # Building the nodes
    # -------------------------------------------------------------------------------------------
    
    # external and internal raidus of the ring
    Rin = 9
    Rout = 10
    Rmid = 0.5*(Rin + Rout)
    nelem = 24
    alpha_div = 2*pi/(nelem)
    
    # positions
    X₀ =  Vec3{T}[]
    
    push!(X₀, Vec3(Rmid, 0, 0))
    
    for idiv = 1:nelem-1
        alphai = idiv*alpha_div
        Xi=Rmid*cos(alphai)
        Zi=Rmid*sin(alphai)
        push!(X₀, Vec3(Xi,  0, Zi))
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
    plane = fill("xz", length(X₀))
    
    # nodes StructArray
    nodes = nodes(X₀, u⁰, u̇⁰, ü⁰, w⁰, ẇ⁰, ẅ⁰, nothing, T)
    
    # -------------------------------------------------------------------------------------------
    # Building the beams
    # -------------------------------------------------------------------------------------------
    
    # total number of beams
    nbeams = nnodes # the ring is closed
    
    # conn
    conn = Vec2{Int}[]
    aux1 = 1:nbeams-1
    aux2 = 2:nbeams
    
    for i in 1:nbeams-1
        push!(conn, (aux1[i], aux2[i])) 
    end
    
    push!(conn, (conn[1][1], conn[end][2]))
    
    # interpolation points per beam
    nbInterpolationPoints = 30
    
    # geometric dimension of the cylinder
    r = 0.5
    E = 100
    ν = 0.0001
    ρ = 0.01
    G = E/(2*(1+ν))
    A = pi*r^2
    I₂₂ = 0.25*pi*r^4
    I₃₃ = 0.25*pi*r^4
    Iₒ = I₂₂+I₃₃
    Iᵣᵣ = Iₒ
    J = Iₒ
    Jᵨ = Mat33(ρ*Iₒ, 0, 0, 0, ρ*I₂₂, 0, 0, 0, ρ*I₃₃)
    Aᵨ = ρ*A
    
    geometry = Geometry{T}(A, I₂₂, I₃₃, Iₒ, Iᵣᵣ, J)
    material = Material{T}(E, G, Aᵨ, Jᵨ)
    
    # beams vector
    beams = beams(nodes, conn, material, geometry, nbInterpolationPoints, nothing)
    
    #-----------------------------------------------------------------------------------
    # Simulation parameters
    # -------------------------------------------------------------------------------------------
    
    # integration parameters
    α = -0.05
    β = 0.25*(1-α)^2
    γ = 0.5*(1-2*α)
    damping = 0
    
    # time step and total time
    Δt = 0.2/16
    Δt_plot = 0.2/16
    tᵉⁿᵈ =  1
    
    # tolerance and maximum number of iterations
    tol_res = 1e-5
    tol_ΔD = 1e-5
    max_it = 10
    
    # Gauss points
    nᴳ = 3
    ωᴳ = Vec3(5/9, 8/9, 5/9)
    zᴳ = Vec3(-sqrt(3/5), 0, sqrt(3/5)) 
    
    # penalty parameters
    kₙ = 5000
    μ = 0
    εᵗ = 0.1
    
    comp = SimulationParameters(α, β, γ, damping,  Δt, Δt_plot, tᵉⁿᵈ, tol_res, tol_ΔD, max_it, nᴳ, ωᴳ, zᴳ, kₙ, μ, εᵗ, T)
    
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
    constraints = nothing
    
    # number of dof (6 per node)
    ndofs = nnodes*6
    
    # Dirichlet boundary conditions: fixed positions
    aux1 = 1:6:ndofs-5
    aux2 = 3:6:ndofs-3
    
    # Fixed dof (x6): I am fixing the r (changes according to the bcs) and z coordinates (blocked)
    fixed_dofs = sort(unique([aux1; aux2])) 
    free_dofs = setdiff(1:ndofs, fixed_dofs) 
    
    # Dirichlet boundary conditions: moving positions
    flag_cylindrical = 1
    disp_dofs = 1:2:(length(fixed_dofs)-1) 
    dispA = -Rmid/2
    tA = 3 
    k = dispA/tA
    Fdisp(t) = (k*t).*(t<=tA)
    flag_disp_vector = false
    disp_vals = T[]
    
    # boundary conditions strucutre 
    bcs = BoundaryConditions(fixed_dofs, free_dofs, flag_cylindrical, flag_disp_vector, Fdisp, disp_vals, disp_dofs, T)
    
    # -------------------------------------------------------------------------------------------
    # SDF
    # -------------------------------------------------------------------------------------------
    
    #there is no contact
    sdf = nothing
    # -------------------------------------------------------------------------------------------
    # Configuration
    # -------------------------------------------------------------------------------------------
    
    # configuration: mesh, external forces and boundary conditions
    conf = Configuration(material, geometry, nnodes, ndofs,  ext_forces, bcs, T)
    
    # -------------------------------------------------------------------------------------------
    # Solve
    # -------------------------------------------------------------------------------------------
    
    params = ParamsTest()
    solver!(nodes, beams, conf, comp, sdf, constraints, params, T)       
    
    # -------------------------------------------------------------------------------------------
    # Test
    # -------------------------------------------------------------------------------------------
    
    pos_matlab = Vec3(7.62779550947649, -7.12963656341915e-10, -2.04386164679398)
    @test isapprox(nodes.X₀[end] + nodes.u[end], pos_matlab; atol=1e-6)
    
end 

test_ring_crimping()

