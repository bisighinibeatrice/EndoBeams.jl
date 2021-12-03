function test_sphere()
    
    # -------------------------------------------------------------------------------------------
    # Building the nodes
    # -------------------------------------------------------------------------------------------
    
    # positions
    L = 10/100; dx = L/10
    
    pos =  Vector{Vec3{T}}()
    
    for x in 0:dx:L
        push!(pos, Vec3(0, 0, x))
    end
    
    # total number of nodes
    nnodes = size(pos,1)
    
    # initial conditions
    u_0 = zeros(size(pos,1)*3)
    udt_0 = zeros(size(pos,1)*3)
    udtdt_0 = zeros(size(pos,1)*3)
    w_0 = zeros(size(pos,1)*3)
    wdt_0 = zeros(size(pos,1)*3)
    wdtdt_0 = zeros(size(pos,1)*3)
    
    # plane for cylindrical coordinates
    plane = fill("xy", length(pos))
    
    # nodes StructArray
    allnodes = constructor_nodes(pos, u_0, udt_0, udtdt_0, w_0, wdt_0, wdtdt_0, plane, Vector{Mat33{T}}(), T)
    
    # -------------------------------------------------------------------------------------------
    # Building the beams
    # -------------------------------------------------------------------------------------------
    
    # total number of beams
    nbeams = nnodes #!!!
    
    # conn
    conn = Vector{Vec2{Int}}()
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
    rho = 7850
    G = E/(2*(1+nu))
    A = pi*r^2
    I22 = 0.25*pi*r^4
    I33 = 0.25*pi*r^4
    Io = I22+I33
    Irr = Io
    J = Io
    Jrho = Mat33{T}(rho*Io, 0, 0, 0, rho*I22, 0, 0, 0, rho*I33)
    Arho = rho*A
    
    geom = Geometry{T}(A, I22, I33, Io, Irr, J)
    mat = Material{T}(E, G, Arho, Jrho)
    
    # beams vector
    allbeams = constructor_beams(allnodes, conn, mat, geom, nbInterpolationPoints, Vector{Mat33{T}}(), T)
    
    #-----------------------------------------------------------------------------------
    # Simulation parameters
    # -------------------------------------------------------------------------------------------
    
    # integration parameters
    alpha = -0.05
    beta = 0.25*(1-alpha)^2
    gamma = 0.5*(1-2*alpha)
    damping = 0
    
    # time step and total time
    dt = 0.01
    dt_plot = 0.1
    tend =  2.5
    
    # tolerance and maximum number of iterations
    tol_res = 1e-5
    tol_ddk = 1e-5
    max_it = 10
    
    # Gaussian points
    nG = 3
    wG = Vec3(5/9, 8/9, 5/9)
    zG = Vec3(-sqrt(3/5), 0, sqrt(3/5)) 
    
    # penalty parameters
    eps_C = 50
    mu_T = 0.3
    eps_tol_fric = 0.1
    
    comp = constructor_simulation_parameters(alpha, beta, gamma, damping,  dt, dt_plot, tend, tol_res, tol_ddk, max_it, nG, wG, zG, eps_C, mu_T, eps_tol_fric, T)
    
    # -------------------------------------------------------------------------------------------
    # External forces
    # -------------------------------------------------------------------------------------------
    
    # external force and applied dof
    flag_crimping = false
    Fext(t) = 0
    dofs_load = []
    
    ext_forces = constructor_ext_force(flag_crimping, Fext, dofs_load, T)
    
    # -------------------------------------------------------------------------------------------
    # Boundary conditions
    # -------------------------------------------------------------------------------------------
    
    # multifreedom constrains
    cons = []
    
    # number of dof (6 per node)
    ndofs = nnodes*6
    
    # Dirichlet boundary conditions: fixed positions
    fixed_dofs = 1:6
    free_dofs = setdiff(1:ndofs, fixed_dofs) 
    
    # Dirichlet boundary conditions: moving positions
    flag_cylindrical = false
    dofs_disp = [1]
    dispA = 10/100
    tA = 2.5
    k = dispA/tA
    Fdisp(t) = (k*t).*(t.<=tA) .+ (-k*t .+ 2*dispA).*((t.>tA) .& (t.<=2*tA)) .+ (k*(t.-2*tA)).*((t.>2*tA) .& (t.<=3*tA)) .+ (-k*t .+ 4*dispA).*((t.>3*tA) .& (t.<=4*tA))
    flag_disp_vector = false
    udisp = []
    
    # boundary conditions strucutre 
    bc = constructor_boundary_conditions(fixed_dofs, free_dofs, flag_cylindrical, flag_disp_vector, Fdisp, udisp, dofs_disp, T)
    
    # -------------------------------------------------------------------------------------------
    # SDF
    # -------------------------------------------------------------------------------------------
    
    # sphere center
    x0 = 0.05
    y0 = 0.005
    z0 = 0.08
    
    # sphere radius
    R = 0.04
    
    sdf = SDF_Sphere{T}(r, R, x0, y0, z0)
    
    # -------------------------------------------------------------------------------------------
    # Configuration
    # -------------------------------------------------------------------------------------------
    
    # configuration: mesh, external forces and boundary conditions
    conf = constructor_configuration(mat, geom, nnodes, ndofs, ext_forces, bc, T)
    
    # -------------------------------------------------------------------------------------------
    # Solve
    # -------------------------------------------------------------------------------------------
    
    params = ParamsTest()
    solver!(allnodes, allbeams, conf, comp, sdf, cons, params, T)         
    
    # -------------------------------------------------------------------------------------------
    # Test 
    # -------------------------------------------------------------------------------------------
    
    rtol = 1e-2
    
    # pos_matlab = Vec3(0.107853165032904, 0.0132224390067658, 0.0360871249073232)
    # t1 = isapprox(allnodes.pos[5] + allnodes.u[5], pos_matlab; rtol)
    # @test t1
    
    # pos_matlab = Vec3(0.129725721716302, 0.0548254754528272, 0.0730574700480992)
    # t2 = isapprox(allnodes.pos[end] + allnodes.u[end], pos_matlab; rtol)
    # @test t2
    
end 

test_sphere()