function test1_ring_plane()
    
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
    pos =  Vec3{T}[]
    push!(pos, Vec3(Rmid, 0, 0))
    
    for idiv = 1:nelem-1
        alphai = idiv*alpha_div
        Xi=Rmid*cos(alphai)
        Zi=Rmid*sin(alphai)
        push!(pos, Vec3(Xi,  0, Zi))
    end 
    
    # total number of nodes
    nnodes = size(pos,1)
    
    # initial conditions
    u_0 = zeros(size(pos,1)*3)
    
    vx_ini = 2*sqrt(2)/2
    vz_ini = -2*sqrt(2)/2
    udt_0 = zeros(size(pos,1)*3,1)
    udt_0[1:3:end-2] .= vx_ini
    udt_0[3:3:end] .= vz_ini
    
    udtdt_0 = zeros(size(pos,1)*3)
    w_0 = zeros(size(pos,1)*3)
    wdt_0 = zeros(size(pos,1)*3)
    wdtdt_0 = zeros(size(pos,1)*3)
    
    # plane for cylindrical coordinates
    plane = fill("xy", length(pos))
    
    # nodes StructArray
    allnodes = constructor_nodes(pos, u_0, udt_0, udtdt_0, w_0, wdt_0, wdtdt_0, plane)
    
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
    nu = 0.0001
    rho = 0.01
    G = E/(2*(1+nu))
    A = pi*r^2
    I22 = 0.25*pi*r^4
    I33 = 0.25*pi*r^4
    Io = I22+I33
    Irr = Io
    J = Io
    Jrho = Mat33(rho*Io, 0, 0, 0, rho*I22, 0, 0, 0, rho*I33)
    Arho = rho*A
    
    geom = Geometry{T}(A, I22, I33, Io, Irr, J)
    mat = Material{T}(E, G, Arho, Jrho)
    
    # beams vector
    allbeams = constructor_beams(allnodes, conn, mat, geom, nbInterpolationPoints)
    
    #-----------------------------------------------------------------------------------
    # Simulation parameters
    # -------------------------------------------------------------------------------------------
    
    # integration parameters
    alpha = -0.05
    beta = 0.25*(1-alpha)^2
    gamma = 0.5*(1-2*alpha)
    damping = 0
    
    # time step and total time
    dt = 0.2/16
    dt_plot = 0.2/16
    tend =  2
    
    # tolerance and maximum number of iterations
    tol_res = 1e-5
    tol_ddk = 1e-5
    max_it = 10
    
    # Gaussian points
    nG = 3
    wG = Vec3(5/9, 8/9, 5/9)
    zG = Vec3(-sqrt(3/5), 0, sqrt(3/5)) 
    
    # penalty parameters
    eps_C = 5000
    mu_T = 0
    eps_tol_fric = 0.1
    
    comp = constructor_simulation_parameters(alpha, beta, gamma, damping,  dt, dt_plot, tend, tol_res, tol_ddk, max_it, nG, wG, zG, eps_C, mu_T, eps_tol_fric)
    
    # -------------------------------------------------------------------------------------------
    # External forces
    # -------------------------------------------------------------------------------------------
    
    # external force and applied dof
    flag_crimping = false
    Fext(t) = 0
    dofs_load = T[]
    
    ext_forces = constructor_ext_force(flag_crimping, Fext, dofs_load)
    
    # -------------------------------------------------------------------------------------------
    # Boundary conditions
    # -------------------------------------------------------------------------------------------
    
    # multifreedom constrains
    cons = T[]
    
    # number of dof (6 per node)
    ndofs = nnodes*6
    
    # Dirichlet boundary conditions: fixed positions
    fixed_dofs = T[]
    free_dofs = setdiff(1:ndofs, fixed_dofs) 
    
    # Dirichlet boundary conditions: moving positions
    flag_cylindrical = false
    dofs_disp = T[]
    Fdisp(t) = 0
    flag_disp_vector = false
    udisp = T[]
    
    # boundary conditions strucutre 
    bc = constructor_boundary_conditions(fixed_dofs, free_dofs, flag_cylindrical, flag_disp_vector, Fdisp, udisp, dofs_disp)
    
    # -------------------------------------------------------------------------------------------
    # SDF
    # -------------------------------------------------------------------------------------------
    
    r = 0.5
    z0 = -10.1
    sdf = SDF_Plane_z{T}(r, z0)
    
    # -------------------------------------------------------------------------------------------
    # Configuration
    # -------------------------------------------------------------------------------------------
    
    # configuration: mesh, external forces and boundary conditions
    conf = constructor_configuration(mat, geom, nnodes, ndofs, ext_forces, bc)
    
    # -------------------------------------------------------------------------------------------
    # Solve
    # -------------------------------------------------------------------------------------------
    
    params = ParamsTest()
    solver!(allnodes, allbeams, conf, comp, sdf, cons, params)       
    
    # -------------------------------------------------------------------------------------------
    # Test 
    # -------------------------------------------------------------------------------------------
    
    pos_matlab = Vec3(13.1644678909684, -1.77821401275396e-14, -4.35229071287855)
    @test isapprox(allnodes.pos[end] + allnodes.u[end], pos_matlab)
    
    pos_matlab = Vec3(-7.30264210560686, -1.74918545324088e-15, 0.536134314694951)
    @test isapprox(allnodes.pos[12] + allnodes.u[12], pos_matlab)
    
end 

test1_ring_plane()