#----------------------------------
# STRUCTURES
#----------------------------------

# Variables for the time and space integration and contact linearization
struct SimulationParameters{T}
    
    # integration parameters
    alpha::T
    beta::T 
    gamma::T
    
    # numerical damping
    damping::T
    
    # time step and total time
    dt::T
    dt_plot::T
    tend::T
    
    # tolerance and maximum number of iterations
    tol_res::T
    tol_ΔDk::T
    max_it::T
    
    # Gaussian points
    nG::Int
    wG::Vec3{T}
    zG::Vec3{T}
    
    # penalty parameters
    epsC::T
    mu_T::T
    eps_tol_fric::T
    
end

# Dirichlet boundary conditions
mutable struct BoundaryConditions{T, F1}
    
    # blocked dofs 
    fixed_dofs::Vector{Int}
    
    # free dofs 
    free_dofs::Vector{Int}
    
    # cylindrical or carthesian coordinate system
    flag_cylindrical::Bool

    # displacement expressed as vector or as function
    flag_disp_vector::Bool
    
    # displacement function
    Fdisp::F1

    # displacement vector
    udisp::Vector{T}

    # dofs where the displacement function is applied
    dof_disps::Vector{Int}   

end

# External forces
struct ExternalForces{T, F}

    # flag crimping: if true, need to convert the force is expressed in radial coordinates and needs to be converted
    flag_crimping::Bool
    
    # concentrated force function 
    Fext::F
    
    # dofs where the concentrated force is applied
    dof_load::Vector{Int}
    
end 

# Material properties
struct Material{T}
    
    E::T
    G::T
    Arho::T
    Jrho::Mat33{T}
    
end 

# Geometrical properties
struct Geometry{T}
    
    A::T
    I22::T
    I33::T
    Io::T
    Irr::T
    J::T   
    
end 

# Configuration of the mesh
struct Configuration{T, F1, F2}
    
    # material
    mat::Material{T}
    
    # geometry
    geom::Geometry{T}
    
    # dof
    ndofs::Int
    disp_dofs::Vector{Int}
    ang_dofs::Vector{Int}
    
    # external forces
    ext_forces::ExternalForces{T,F1}
    
    # boundary conditions
    bc::BoundaryConditions{T,F2}
    
end

# Force vectors needed at the next time step by the solver
struct Solution{T}  
    
    fext::Vector{T}
    Tint::Vector{T}
    Tk::Vector{T}
    Tct::Vector{T}
    Tconstr::Vector{T}
    
end 

# Current nodal solutions (preallocation)
struct NodalSolution{T}  
    
    D::Vector{T}
    Ddt::Vector{T}
    Ddtdt::Vector{T}
    Ddt_n::Vector{T}
    
    r::Vector{T}
    aux1::Vector{T} #predictor
    aux2::Vector{T} #predictor
    Ktan::SparseMatrixCSC{T,Int}
    ΔD::Vector{T}    
    asol::Vector{T}
    f_aux::Vector{T}

    r_free::Vector{T}
    Ktan_free::SparseMatrixCSC{T,Int}
    ΔD_free::Vector{T}
    sparsity_map_free::Vector{Int}  
    
end 

# Global matrices structure (sparse arrays)
struct Matrices{T}
    
    Kint::SparseMatrixCSC{T,Int}
    Ck::SparseMatrixCSC{T,Int} 
    M::SparseMatrixCSC{T,Int}
    Kct::SparseMatrixCSC{T,Int}
    Kconstr::SparseMatrixCSC{T,Int}
    Cconstr::SparseMatrixCSC{T,Int}
    
    Tint::Vector{T}
    Tk::Vector{T}
    Tdamp::Vector{T}
    Tct::Vector{T}
    Tconstr::Vector{T}
    
end 

# Energy contributions structure
mutable struct Energy{T}
    
    Phi_energy::T
    K_energy::T
    C_energy::T
    
end

# Fixed matrices used in the beam contributions computation (preallocation)
struct PreAllocatedMatricesFixed{T}
    
    E1::Vec3{T}
    E2::Vec3{T}
    E3::Vec3{T}
    
    A1::Mat312{T}
    A2::Mat312{T}
    re::Vec12{T}
    
    P1G_v::Vec3{Mat36{T}}
    P2G_v::Vec3{Mat36{T}}
    NG_v::Vec3{Mat312{T}}
    
end 

# Information related to the contact at the Gaussian points
mutable struct GPSolution{T}  
    
    status::Vector{Int}
    xGP::Vector{Vec3{T}}
    fGP_N::Vector{Vec3{T}}
    fGP_T::Vector{Vec3{T}}
    gGP::Vector{T}
    
end

#----------------------------------
# PUBLIC CONSTRUCTORS
#----------------------------------

"""
Constructor of the structure containing the material properties

Material{T}(E, G, Arho, Jrho) =
constructor_material_properties(E, nu, rho, rWireSection)

Returns a Material structure.

"""
function constructor_material_properties(E, nu, rho, rWireSection)

    G = E/(2*(1+nu))
    A = pi*rWireSection^2
    I22 = 0.25*pi*rWireSection^4
    I33 = 0.25*pi*rWireSection^4
    Io = I22+I33
    Jrho = Mat33{T}(rho*Io, 0, 0, 0, rho*I22, 0, 0, 0, rho*I33)
    Arho = rho*A
 
    mat = Material{T}(E, G, Arho, Jrho)

    return mat

end

"""
Constructor of the structure containing the geometrical properties.

Geometry{T}(A, I22, I33, Io, Irr, J) = 
constructor_geometry_properties(rWireSection)

Returns a Geometry structure.

"""
function constructor_geometry_properties(rWireSection)

    A = pi*rWireSection^2
    I22 = 0.25*pi*rWireSection^4
    I33 = 0.25*pi*rWireSection^4
    Io = I22+I33
    Irr = Io
    J = Io
    
    geom = Geometry{T}(A, I22, I33, Io, Irr, J)

    return geom

end 

"""
Constructor of the structure containing the simulation parameters

SimulationParameters{T}(alpha, beta, gamma, damping, dt, dt_plot, tend, tol_res, tol_ddk, max_it, nG, wG, zG, eps_C, mu_T, eps_tol_fric) = 
constructor_simulation_parameters(alpha, beta, gamma, damping, dt, dt_plot, tend, tol_res, tol_ddk, max_it, nG, wG, zG, eps_C, mu_T, eps_tol_fric, T=Float64)    

Returns a SimulationParameters structure.

"""

function constructor_simulation_parameters(alpha, beta, gamma, damping, dt, dt_plot, tend, tol_res, tol_ddk, max_it, nG, wG, zG, eps_C, mu_T, eps_tol_fric, T=Float64)
    
    return SimulationParameters{T}(alpha, beta, gamma, damping, dt, dt_plot, tend, tol_res, tol_ddk, max_it, nG, wG, zG, eps_C, mu_T, eps_tol_fric)
    
end 

"""
Constructor of the structure containing the information about the external load, if present

ExternalForces{T, typeof(Fext)}(flag_crimping, Fext, dof_load) = 
constructor_ext_force(flag_crimping, Fext, dof_load, T=Float64)

Returns a ExternalForces structure.

"""
function constructor_ext_force(flag_crimping, Fext, dof_load, T=Float64)
        
    return ExternalForces{T, typeof(Fext)}(flag_crimping, Fext, dof_load)
    
end 

"""
Constructor of the structure containing the information about the boundary conditions, if present

BoundaryConditions{T, F1}(fixed_dofs, free_dofs, flag_cylindrical, flag_disp_vector, Fdisp, udisp, dofs_disp) = 
constructor_boundary_conditions(fixed_dofs, free_dofs, flag_cylindrical, flag_disp_vector, Fdisp, udisp, dofs_disp, T=Float64)

Returns a BoundaryConditions structure.

"""
function constructor_boundary_conditions(fixed_dofs, free_dofs, flag_cylindrical, flag_disp_vector, Fdisp, udisp, dofs_disp, T=Float64)
    
    F1 = typeof(Fdisp)

    return  BoundaryConditions{T, F1}(fixed_dofs, free_dofs, flag_cylindrical, flag_disp_vector, Fdisp, udisp, dofs_disp)
    
end 

"""
Configuration{T, typeof(ext_forces.Fext), typeof(bc.Fdisp)}(mat, geom, ndofs, disp_dof, ang_dof, ext_forces, bc) = 
constructor_configuration(mat, geom, nnodes, ndofs, ext_forces, bc, T=Float64)

Returns a Configuration structure.

"""
function constructor_configuration(mat, geom, nnodes, ndofs, ext_forces, bc, T=Float64)
    
    # displacement dofs
    disp_dofs = zeros(3*nnodes)
    disp_dofs[1:3:end-2] = 1:6:ndofs-5
    disp_dofs[2:3:end-1] = 2:6:ndofs-4
    disp_dofs[3:3:end] = 3:6:ndofs-3

    # rotation dofs
    ang_dofs = setdiff(1:ndofs,disp_dofs)

    return Configuration{T, typeof(ext_forces.Fext), typeof(bc.Fdisp)}(mat,geom,ndofs, disp_dofs, ang_dofs, ext_forces,  bc)
    
end 

#----------------------------------
# PRIVATE CONSTRUCTORS
#----------------------------------

# Constructor of the structure where the last force vectors are saved in order to be used in the next step by the solver
function constructor_solution(conf, T=Float64)
    
    ndofs = conf.ndofs

    # external force @t=0
    fext_0  = zeros(ndofs)
    get_current_external_force!(fext_0, 0, conf)
    
    Tint_0 = zeros(ndofs)
    Tk_0 = zeros(ndofs)
    Tct_0 = zeros(ndofs)
    Tconstr_0 = zeros(ndofs)
    
    return Solution{T}(fext_0, Tint_0, Tk_0, Tct_0, Tconstr_0)
    
end 

# Constructor of the structure containing some fixed matrices used in the beam contributions computation 
function constructor_preallocated_matrices_fixed(allbeams, comp, T=Float64)
    
    E1 = Vec3(1,0,0)
    E2 = Vec3(0,1,0)
    E3 = Vec3(0,0,1)
    
    # Note: do not use these matrices, define a function that reorganises elements insead
    A1 = Mat312(
    0, 0, 0, 
    0,-1, 0, 
    0, 0,-1, 
    0, 0, 0,
    0, 0, 0,
    0, 0, 0,
    0, 0, 0,
    0, 1, 0,
    0, 0, 1,
    0, 0, 0,
    0, 0, 0,
    0, 0, 0)
    
    A2 = Mat312(
    0, 0, 0, 
    0, 0,-1, 
    0, 1, 0, 
    0, 0, 0,
    0, 0, 0,
    0, 0, 0,
    0, 0, 0,
    0, 0, 1,
    0,-1, 0,
    0, 0, 0,
    0, 0, 0,
    0, 0, 0)
    
    re = Vec12( 
    -1,
    0,
    0,
    0,
    0,
    0,
    1,
    0,
    0,
    0,
    0,
    0 )
    
    # ----------------------------------------------------------------------
    l0 = allbeams.l0[1]
    
    # _1
    zG = comp.zG[1] 
    xG = l0*(zG+1)/2
    
    # eqB3:5 in [2]: shape functions
    N1_1 = 1-xG/l0
    N2_1 = 1-N1_1
    N3_1 = xG*(1-xG/l0)^2
    N4_1 = -(1-xG/l0)*((xG^2)/l0)
    N5_1 = (1-3*xG/l0)*(1-xG/l0)
    N6_1 = (3*xG/l0-2)*(xG/l0)
    N7_1 = N3_1+N4_1
    N8_1 = N5_1+N6_1-1
    
    # eqB1 in [2]: FE interpolation matrices for displacement
    P1G_1 = Mat36(
    0, 0, 0, 
    0, 0, -N3_1, 
    0, N3_1, 0, 
    0, 0, 0,
    0, 0, -N4_1,
    0, N4_1, 0)
    
    # eqB2 in [2]: FE interpolation matrices for rotations
    P2G_1 = Mat36(
    N1_1, 0, 0, 
    0, N5_1, 0, 
    0, 0, N5_1, 
    N2_1, 0, 0,
    0, N6_1, 0,
    0, 0, N6_1)
    
    # eqD7 in [2], N
    NG_1 = Mat312(
    N1_1, 0, 0, 
    0, N1_1, 0, 
    0, 0, N1_1, 
    0, 0, 0,
    0, 0, 0,
    0, 0, 0,
    N2_1, 0, 0,
    0, N2_1, 0,
    0, 0, N2_1,
    0, 0, 0,
    0, 0, 0,
    0, 0, 0)
    
    # NG = Vec4{T}(N1*eye33, zeros(Mat33{T}), N2*eye33, zeros(Mat33{T}))
    
    # ----------------------------------------------------------------------
    # _2
    zG = comp.zG[2] 
    xG = l0*(zG+1)/2
    
    # eqB3:5 in [2]: shape functions
    N1_2 = 1-xG/l0
    N2_2 = 1-N1_2
    N3_2 = xG*(1-xG/l0)^2
    N4_2 = -(1-xG/l0)*((xG^2)/l0)
    N5_2 = (1-3*xG/l0)*(1-xG/l0)
    N6_2 = (3*xG/l0-2)*(xG/l0)
    N7_2 = N3_2+N4_2
    N8_2 = N5_2+N6_2-1
    
    # eqB1 in [2]: FE interpolation matrices for displacement
    P1G_2 = Mat36(
    0, 0, 0, 
    0, 0, -N3_2, 
    0, N3_2, 0, 
    0, 0, 0,
    0, 0, -N4_2,
    0, N4_2, 0)
    
    # eqB2 in [2]: FE interpolation matrices for rotations
    P2G_2 = Mat36(
    N1_2, 0, 0, 
    0, N5_2, 0, 
    0, 0, N5_2, 
    N2_2, 0, 0,
    0, N6_2, 0,
    0, 0, N6_2)
    
    # eqD7 in [2], N
    NG_2 = Mat312(
    N1_2, 0, 0, 
    0, N1_2, 0, 
    0, 0, N1_2, 
    0, 0, 0,
    0, 0, 0,
    0, 0, 0,
    N2_2, 0, 0,
    0, N2_2, 0,
    0, 0, N2_2,
    0, 0, 0,
    0, 0, 0,
    0, 0, 0)
    
    # ----------------------------------------------------------------------
    # _3
    zG = comp.zG[3] 
    xG = l0*(zG+1)/2
    
    # eqB3:5 in [2]: shape functions
    N1_3 = 1-xG/l0
    N2_3 = 1-N1_3
    N3_3 = xG*(1-xG/l0)^2
    N4_3 = -(1-xG/l0)*((xG^2)/l0)
    N5_3 = (1-3*xG/l0)*(1-xG/l0)
    N6_3 = (3*xG/l0-2)*(xG/l0)
    N7_3 = N3_3+N4_3
    N8_3 = N5_3+N6_3-1
    
    # eqB1 in [2]: FE interpolation matrices for displacement
    P1G_3 = Mat36(
    0, 0, 0, 
    0, 0, -N3_3, 
    0, N3_3, 0, 
    0, 0, 0,
    0, 0, -N4_3,
    0, N4_3, 0)
    
    
    # eqB2 in [2]: FE interpolation matrices for rotations
    P2G_3 = Mat36(
    N1_3, 0, 0, 
    0, N5_3, 0, 
    0, 0, N5_3, 
    N2_3, 0, 0,
    0, N6_3, 0,
    0, 0, N6_3)
    
    
    # eqD7 in [2], N
    NG_3 = Mat312(
    N1_3, 0, 0, 
    0, N1_3, 0, 
    0, 0, N1_3, 
    0, 0, 0,
    0, 0, 0,
    0, 0, 0,
    N2_3, 0, 0,
    0, N2_3, 0,
    0, 0, N2_3,
    0, 0, 0,
    0, 0, 0,
    0, 0, 0)
    
    # ----------------------------------------------------------------------
    # Assembly
    
    P1G_v = Vec3{Mat36{T}}(P1G_1, P1G_2, P1G_3)
    P2G_v = Vec3{Mat36{T}}(P2G_1, P2G_2, P2G_3)
    NG_v = Vec3{Mat312{T}}(NG_1, NG_2, NG_3)
    return PreAllocatedMatricesFixed{T}(E1, E2, E3, A1, A2, re, P1G_v, P2G_v, NG_v)
    
end 

# Constructor of the structure containing the energy contributions
function constructor_energy(T=Float64)
    
    return Energy{T}(0, 0, 0)
    
end

# Constructor of the structure containing the preallocated variables used in the solver
function constructor_nodal_solution(ndofs, nfreedofs, dofs_tan, dofs_free, spmap_free, T=Float64)
    
    rows_tan, cols_tan, zval_tan = dofs_tan
    rows_free, cols_free, zval_free = dofs_free

    return NodalSolution{T}(

    zeros(ndofs), #D
    zeros(ndofs), #Ddt 
    zeros(ndofs), #Ddtdt 
    zeros(ndofs), #Ddt_n  -> needed by predictor
    
    spzeros(ndofs), #r
    spzeros(ndofs), #res -> needed by predictor
    spzeros(ndofs), #res -> needed by predictor
    sparse(rows_tan, cols_tan, zval_tan), #Ktan
    zeros(ndofs), #ΔD    
    zeros(ndofs), #asol
    zeros(ndofs), #f_aux

    zeros(nfreedofs), #r_free
    sparse(rows_free, cols_free, zval_free), #Ktan_free
    zeros(nfreedofs), #ΔD_free
    spmap_free)

end

# Constructor of the structure containing the global matrices 
function constructor_global_matrices(allnodes, dofs_tan, T=Float64)
    
    rows, cols, zval = dofs_tan

    nnodes = length(allnodes)
    ndofs_per_node = 6
    ndofs = nnodes*ndofs_per_node
    
    Kint_0 = sparse(rows, cols, zval)
    Ck_0 =  sparse(rows, cols, zval)
    M_0 = sparse(rows, cols, zval)
    Kct_0 =  sparse(rows, cols, zval)
    Kconstr_0  = sparse(rows, cols, zval)
    Ccosntr_0  = sparse(rows, cols, zval)

    Tint_0 = zeros(ndofs)
    Tk_0 =  zeros(ndofs)
    Tdamp_0 = zeros(ndofs)
    Tct_0 = zeros(ndofs)
    Tconstr_0  = spzeros(ndofs)

    return Matrices{T}(Kint_0, Ck_0, M_0, Kct_0, Kconstr_0, Ccosntr_0,Tint_0, Tk_0, Tdamp_0, Tct_0, Tconstr_0)
    
end 

# Constructor of the sparse matrices
function constructor_sparse_matrices!(allbeams, allnodes, pncons, conf, T=Float64)

    dofs_tan, dofs_free, spmap_free = compute_sparsity!(allbeams, allnodes, pncons, conf)
    
    ndofs = length(allnodes)*6
    nfreedofs = length(setdiff(1:ndofs, conf.bc.fixed_dofs))

    matrices = constructor_global_matrices(allnodes, dofs_tan)
    nodes_sol = constructor_nodal_solution(ndofs, nfreedofs, dofs_tan, dofs_free, spmap_free)

    return matrices, nodes_sol

end 

# Constructor of the structure where the information about the GP are saved
function constructor_solution_GP(nbeams, T=Float64)
    
    nGP = nbeams*3
    xGP = Vector{Vec3{T}}()
    fGP_N = Vector{Vec3{T}}()
    tGP_T = Vector{Vec3{T}}()
    for i in 1:nGP
        push!(xGP, zeros(Vec3{T}))
        push!(fGP_N, zeros(Vec3{T}))
        push!(tGP_T, zeros(Vec3{T}))
    end
    gGP = zeros(nGP)
    status = zeros(nGP)
    
    return GPSolution{T}(status, xGP, fGP_N, tGP_T, gGP)
    
end 