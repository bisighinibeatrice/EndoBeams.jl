#----------------------------------
# STRUCTURES
#----------------------------------

# Variables for the time and space integration and contact linearization
struct SimulationParameters{T}
    
    # integration parameters
    α::T 
    β::T 
    γ::T 
    
    # numerical damping
    damping::T
    
    # time step and total time
    Δt::T
    Δt_plot::T
    tᵉⁿᵈ::T
    
    # tolerance and maximum number of iterations
    res_tol::T
    tol_ΔDk::T
    max_it::T
    
    # Gauss points
    nᴳ::Int
    ωᴳ::Vec3{T}
    zᴳ::Vec3{T}
    
    # penalty parameters
    εᶜ::T
    μ::T
    εᵗ::T
    
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
    Aᵨ::T
    Jᵨ::Mat33{T}
    
end 

# Geometrical properties
struct Geometry{T}
    
    A::T
    I₂₂::T
    I₃₃::T
    Iₒ::T
    Iᵣᵣ::T
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
    
    fᵉˣᵗ::Vector{T}
    Tⁱⁿᵗ::Vector{T}
    Tᵏ::Vector{T}
    Tᶜ::Vector{T}
    Tᶜᵒⁿ::Vector{T}
    
end 

# Current nodal solutions (preallocation)
struct NodalSolution{T}  
    
    D::Vector{T}
    Ḋ::Vector{T}
    D̈::Vector{T}
    
    r::Vector{T}
    Ktan_mat::SparseMatrixCSC{T,Int}
    Ktan::Vector{Float64}
    ΔD::Vector{T}    
    temp::Vector{T}

    r_free::Vector{Float64}
    Ktan_free::SparseMatrixCSC{Float64,Int}
    ΔD_free::Vector{Float64}
    sparsity_map_free::Vector{Int}  
    
end 

# Global matrices structure (sparse arrays)
struct Matrices{T}
    
    K_mat ::SparseMatrixCSC{T,Int}
    K::Vector{T}
    K_mtbufs::Vector{Vector{T}}
    C_mat::SparseMatrixCSC{T,Int}
    C::Vector{T}
    C_mtbufs::Vector{Vector{T}}
    M_mat::SparseMatrixCSC{T,Int}
    M::Vector{T}
    M_mtbufs::Vector{Vector{T}}
    
    Tⁱⁿᵗ::Vector{T}
    Tⁱⁿᵗ_mtbufs::Vector{Vector{T}}
    Tᵏ::Vector{T}
    Tᵏ_mtbufs::Vector{Vector{T}}
    Tᶜ::Vector{T}
    Tᶜ_mtbufs::Vector{Vector{T}}
    Tᶜᵒⁿ::Vector{T}
    
end 

# Energy contributions structure
mutable struct Energy{T}
    
    strain_energy::T
    kinetic_energy::T
    contact_energy::T
    
end


# Information related to the contact at the Gauss points
struct GPSolution{T}  
    
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
    mat = constructor_material_properties(E, nu, ρ, rWireSection)

Constructor of the structure containing the material properties:
- `E`: Young modulus;
- `nu`: Poisson coefficient;
- `ρ`: density;
- `rWireSection`: beam radius.

Returns a Material structure.
"""
function constructor_material_properties(E, nu, ρ, rWireSection, T=Float64)

    G = E/(2*(1+nu))
    A = pi*rWireSection^2
    I₂₂ = pi*rWireSection^4/4
    I₃₃ = pi*rWireSection^4/4
    Iₒ = I₂₂+I₃₃
    Jᵨ = Mat33(ρ*Iₒ, 0, 0, 0, ρ*I₂₂, 0, 0, 0, ρ*I₃₃)
    Aᵨ = ρ*A

    mat = Material{T}(E, G, Aᵨ, Jᵨ)

    return mat

end

"""
    geom = constructor_geometry_properties(rWireSection)

Constructor of the structure containing the geometrical properties:
- `rWireSection`: beam radius.

Returns a Geometry structure.
"""
function constructor_geometry_properties(rWireSection, T=Float64)

    A = pi*rWireSection^2
    I₂₂ = pi*rWireSection^4/4
    I₃₃ = pi*rWireSection^4/4
    Iₒ = I₂₂+I₃₃
    Iᵣᵣ = Iₒ
    J = Iₒ
    
    geom = Geometry{T}(A, I₂₂, I₃₃, Iₒ, Iᵣᵣ, J)

    return geom

end 

"""
    comp = constructor_simulation_parameters(α, β, γ, damping, Δt, Δt_plot, tᵉⁿᵈ, res_tol, tol_ddk, max_it, nG, ωG, zG, εᶜ, μ, εᵗ, T=Float64)    

Constructor of the structure containing the simulation parameters:
    - `α`: integration parameter;
    - `β`: integration parameter;
    - `γ`: integration parameter;
    - `damping`: damping coefficient;
    - `Δt`: time step;
    - `Δt_plot`: time step for the saving of output files;
    - `tᵉⁿᵈ`: total time of the simulation;
    - `res_tol`: residual tolerance for the Newton-Raphson algorithm;
    - `tol_ddk`: solution vector tolerance for the Newton-Raphson algorithm;
    - `max_it`: maximum number of iterations for the Newton-Raphson algorithm;
    - `nG`: number of Gauss points;
    - `ωG`: Gauss points weights;
    - `zG`: Gauss points positions along centreline;
    - `εᶜ`: penalty coefficient for the contact normal contributions;
    - `μ`: friction coefficient;
    - `εᵗ`: regularisation coefficient for the contact tangential contributions.

Returns a SimulationParameters structure.
"""
function constructor_simulation_parameters(α, β, γ, damping, Δt, Δt_plot, tᵉⁿᵈ, res_tol, tol_ddk, max_it, nG, ωG, zG, εᶜ, μ, εᵗ, T=Float64)
    
    return SimulationParameters{T}(α, β, γ, damping, Δt, Δt_plot, tᵉⁿᵈ, res_tol, tol_ddk, max_it, nG, ωG, zG, εᶜ, μ, εᵗ)
    
end 

"""
    fᵉˣᵗ = constructor_ext_force(flag_crimping, Fext, dof_load, T=Float64)

Constructor of the structure containing the information about the external load, if present:
- `flag_crimping`: true if the force is defined in cylindral coordinates;
- `Fext`: external load amplitude;
- `dof_load`: DOFs interested by the external load.

Returns a ExternalForces structure.
"""
function constructor_ext_force(flag_crimping, Fext, dof_load, T=Float64)
        
    return ExternalForces{T, typeof(Fext)}(flag_crimping, Fext, dof_load)
    
end 

"""
    bcs = constructor_boundary_conditions(fixed_dofs, free_dofs, flag_cylindrical, flag_disp_vector, Fdisp, udisp, dofs_disp, T=Float64)

Constructor of the structure containing the information about the Dirichlet boundary conditions (BCs), if present:
- `fixed_dofs`: fixed DOFs;
- `free_dofs`: free DOFs;
- `flag_cylindrical`: true if the BCs are defined in cylindral coordinates;
- `flag_disp_vector`: true if the BCs are defined as a vector;
- `Fdisp`: imposed dispacement amplitude;
- `udisp`: vector; containing the imposed dispacements;
- `dofs_disp`: DOFs interested by the Dirichlet BCs.

Returns a BoundaryConditions structure.
"""
function constructor_boundary_conditions(fixed_dofs, free_dofs, flag_cylindrical, flag_disp_vector, Fdisp, udisp, dofs_disp, T=Float64)
    
    F1 = typeof(Fdisp)

    return  BoundaryConditions{T, F1}(fixed_dofs, free_dofs, flag_cylindrical, flag_disp_vector, Fdisp, udisp, dofs_disp)
    
end 

"""
    conf = constructor_configuration(mat, geom, nnodes, ndofs, ext_forces, bcs, T=Float64)

Constructor of the structure collecting the information for the simulation:
- `mat`: material properties (Material{T});
- `geom`: geoltrical properties (Geometry{T});
- `nnodes`: number of nodes in the system;
- `ndofs`: number of DOFs in the system;
- `bcs`: Dirichlet BCs.

Returns a Configuration structure.
"""
function constructor_configuration(mat, geom, nnodes, ndofs, ext_forces, bcs, T=Float64)
    
    # displacement dofs
    disp_dofs = zeros(3*nnodes)
    disp_dofs[1:3:end-2] = 1:6:ndofs-5
    disp_dofs[2:3:end-1] = 2:6:ndofs-4
    disp_dofs[3:3:end] = 3:6:ndofs-3

    # rotation dofs
    ang_dofs = setdiff(1:ndofs,disp_dofs)

    return Configuration{T, typeof(ext_forces.Fext), typeof(bcs.Fdisp)}(mat,geom,ndofs, disp_dofs, ang_dofs, ext_forces, bcs)
    
end 

#----------------------------------
# PRIVATE CONSTRUCTORS
#----------------------------------

# Constructor of the structure where the last force vectors are saved in order to be used in the next step by the solver
function constructor_solution(conf, T=Float64)
    
    ndofs = conf.ndofs

    # external force @t=0
    fᵉˣᵗ_0  = zeros(ndofs)
    get_current_external_force!(fᵉˣᵗ_0, 0, conf)
    
    T⁰ᵢₙₜ = zeros(ndofs)
    T⁰ₖ = zeros(ndofs)
    T⁰ₜ = zeros(ndofs)
    T⁰ₓ = zeros(ndofs)
    
    return Solution{T}(fᵉˣᵗ_0, T⁰ᵢₙₜ, T⁰ₖ, T⁰ₜ, T⁰ₓ)
    
end 


# Constructor of the structure containing the energy contributions
function constructor_energy(T=Float64)
    
    return Energy{T}(0, 0, 0)
    
end

# Constructor of the structure containing the preallocated variables used in the solver
function constructor_nodal_solution(ndofs, nfreedofs, dofs_tan, dofs_free, spmap_free, T=Float64)
    
    rows_tan, cols_tan, nzval_tan = dofs_tan
    rows_free, cols_free, nzval_free = dofs_free

    Ktan_mat = sparse(rows_tan, cols_tan, nzval_tan)

    return NodalSolution{T}(

    zeros(ndofs), #D
    zeros(ndofs), #Ḋ 
    zeros(ndofs), #D̈ 
    
    zeros(ndofs), #r
    Ktan_mat, #Ktan_mat
    nonzeros(Ktan_mat), #Ktan
    zeros(ndofs), #ΔD    
    zeros(ndofs), #temp vector

    zeros(nfreedofs), #r_free
    sparse(rows_free, cols_free, nzval_free), #Ktan_free
    zeros(nfreedofs), #ΔD_free
    spmap_free)

end

# Constructor of the structure containing the global matrices 
function constructor_global_matrices(nodes, dofs_tan, T=Float64)
    
    rows, cols, nzval = dofs_tan
    ncons = length(nzval)

    nnodes = length(nodes)
    ndofs_per_node = 6
    ndofs = nnodes*ndofs_per_node
    
    C_mat = sparse(rows, cols, zeros(ncons))
    C = nonzeros(C_mat)
    C_mtbufs = [copy(C) for _ in 1:Threads.nthreads()]
    
    M_mat = sparse(rows, cols, zeros(ncons))
    M = nonzeros(M_mat)
    M_mtbufs = [copy(M) for _ in 1:Threads.nthreads()]

    K_mat = sparse(rows, cols, zeros(ncons))
    K =  nonzeros(K_mat)
    K_mtbufs = [copy(K) for _ in 1:Threads.nthreads()]

    Tⁱⁿᵗ = zeros(ndofs)
    Tⁱⁿᵗ_mtbufs = [copy(Tⁱⁿᵗ) for _ in 1:Threads.nthreads()]
    Tᵏ =  zeros(ndofs)
    Tᵏ_mtbufs = [copy(Tᵏ) for _ in 1:Threads.nthreads()]
    Tᶜ = zeros(ndofs)
    Tᶜ_mtbufs = [copy(Tᶜ) for _ in 1:Threads.nthreads()]
    Tᶜᵒⁿ = zeros(ndofs)

    return Matrices{T}(K_mat, K, K_mtbufs, C_mat, C, C_mtbufs, M_mat, M, M_mtbufs, Tⁱⁿᵗ, Tⁱⁿᵗ_mtbufs, Tᵏ, Tᵏ_mtbufs, Tᶜ, Tᶜ_mtbufs, Tᶜᵒⁿ)
    
end 

# Constructor of the sparse matrices
function constructor_sparse_matrices!(beams, nodes, pncons, conf, T=Float64)

    dofs_tan, dofs_free, spmap_free = compute_sparsity!(beams, nodes, pncons, conf, T)
    
    ndofs = length(nodes)*6
    nfreedofs = length(setdiff(1:ndofs, conf.bc.fixed_dofs))

    matrices = constructor_global_matrices(nodes, dofs_tan, T)
    nodes_sol = constructor_nodal_solution(ndofs, nfreedofs, dofs_tan, dofs_free, spmap_free, T)

    return matrices, nodes_sol

end 

# Constructor of the structure where the information about the GP are saved
function constructor_solution_GP(nbeams, T=Float64)
    
    nGP = nbeams*3
    xGP = Vec3{T}[]
    fGP_N = Vec3{T}[]
    tGP_T = Vec3{T}[]
    for i in 1:nGP
        push!(xGP, zeros(Vec3{T}))
        push!(fGP_N, zeros(Vec3{T}))
        push!(tGP_T, zeros(Vec3{T}))
    end
    gGP = zeros(nGP)
    status = zeros(nGP)
    
    return GPSolution{T}(status, xGP, fGP_N, tGP_T, gGP)
    
end 
