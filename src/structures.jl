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
    tol_res::T
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
mutable struct BoundaryConditions{T, TF}
    
    # blocked dofs 
    fixed_dofs::Vector{Int}
    
    # free dofs 
    free_dofs::Vector{Int}
    
    # displacement function
    u::TF

    # displacement vector
    disp_vals::Vector{T}

    # dofs where the displacement function is applied
    disp_dofs::Vector{Int}   

end

# External forces
struct ExternalForces{TF}
    
    # concentrated force function 
    f::TF
    
    # dofs where the concentrated force is applied
    loaded_dofs::Vector{Int}
    
end 

# Material properties
struct Material{T, TJ}
    
    E::T
    G::T
    Aᵨ::T
    Jᵨ::TJ
    
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
struct Configuration{T, TJ, TF, TU}
    
    # material
    mat::Material{T, TJ}
    
    # geometry
    geom::Geometry{T}
    
    # dof
    ndofs::Int
    disp_dofs::Vector{Int}
    rot_dofs::Vector{Int}
    
    # external forces
    ext_forces::ExternalForces{TF}
    
    # boundary conditions
    bcs::BoundaryConditions{T, TU}
    
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
    
end 

# Global matrices structure (sparse arrays)
struct Matrices{T}
    
    K_mat ::SparseMatrixCSC{T,Int}
    K::Vector{T}
    C_mat::SparseMatrixCSC{T,Int}
    C::Vector{T}
    M_mat::SparseMatrixCSC{T,Int}
    M::Vector{T}
    
    Tⁱⁿᵗ::Vector{T}
    Tᵏ::Vector{T}
    Tᶜ::Vector{T}
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
    mat = constructor_material_properties(E, nu, ρ, radius)

Constructor of the structure containing the material properties:
- `E`: Young modulus;
- `nu`: Poisson coefficient;
- `ρ`: density;
- `radius`: beam radius.

Returns a Material structure.
"""
function constructor_material_properties(E, nu, ρ, radius, T=Float64)

    G = E/(2*(1+nu))
    A = pi*radius^2
    I₂₂ = pi*radius^4/4
    I₃₃ = pi*radius^4/4
    Iₒ = I₂₂ + I₃₃
    Jᵨ = ρ * Diagonal(Vec3(Iₒ, I₂₂, I₃₃))
    Aᵨ = ρ*A

    mat = Material{T, typeof(Jᵨ)}(E, G, Aᵨ, Jᵨ)

    return mat

end

"""
    geom = constructor_geometry_properties(radius)

Constructor of the structure containing the geometrical properties:
- `radius`: beam radius.

Returns a Geometry structure.
"""
function constructor_geometry_properties(radius, T=Float64)

    A = pi*radius^2
    I₂₂ = pi*radius^4/4
    I₃₃ = pi*radius^4/4
    Iₒ = I₂₂+I₃₃
    Iᵣᵣ = Iₒ
    J = Iₒ
    
    geom = Geometry{T}(A, I₂₂, I₃₃, Iₒ, Iᵣᵣ, J)

    return geom

end 

"""
    comp = constructor_simulation_parameters(α, β, γ, damping, Δt, Δt_plot, tᵉⁿᵈ, tol_res, tol_ΔD, max_it, nG, ωG, zG, εᶜ, μ, εᵗ, T=Float64)    

Constructor of the structure containing the simulation parameters:
    - `α`: integration parameter;
    - `β`: integration parameter;
    - `γ`: integration parameter;
    - `damping`: damping coefficient;
    - `Δt`: time step;
    - `Δt_plot`: time step for the saving of output files;
    - `tᵉⁿᵈ`: total time of the simulation;
    - `tol_res`: residual tolerance for the Newton-Raphson algorithm;
    - `tol_ΔD`: solution vector tolerance for the Newton-Raphson algorithm;
    - `max_it`: maximum number of iterations for the Newton-Raphson algorithm;
    - `nG`: number of Gauss points;
    - `ωG`: Gauss points weights;
    - `zG`: Gauss points positions along centreline;
    - `εᶜ`: penalty coefficient for the contact normal contributions;
    - `μ`: friction coefficient;
    - `εᵗ`: regularisation coefficient for the contact tangential contributions.

Returns a SimulationParameters structure.
"""
function constructor_simulation_parameters(α, β, γ, damping, Δt, Δt_plot, tᵉⁿᵈ, tol_res, tol_ΔD, max_it, nG, ωG, zG, εᶜ, μ, εᵗ, T=Float64)
    
    return SimulationParameters{T}(α, β, γ, damping, Δt, Δt_plot, tᵉⁿᵈ, tol_res, tol_ΔD, max_it, nG, ωG, zG, εᶜ, μ, εᵗ)
    
end 

"""
    fᵉˣᵗ = constructor_ext_force(F, loaded_dofs, T=Float64)

Constructor of the structure containing the information about the external load, if present:
- `F`: external load amplitude;
- `loaded_dofs`: DOFs interested by the external load.

Returns a ExternalForces structure.
"""
function constructor_ext_force(force_fun::T, loaded_dofs) where T
        
    return ExternalForces{T}(force_fun, loaded_dofs)
    
end 

"""
    bcs = constructor_boundary_conditions(fixed_dofs, free_dofs, flag_cylindrical, flag_disp_vector, Fdisp, disp_vals, disp_dofs, T=Float64)

Constructor of the structure containing the information about the Dirichlet boundary conditions (BCs), if present:
- `fixed_dofs`: fixed DOFs;
- `free_dofs`: free DOFs;
- `flag_cylindrical`: true if the BCs are defined in cylindral coordinates;
- `flag_disp_vector`: true if the BCs are defined as a vector;
- `Fdisp`: imposed dispacement amplitude;
- `disp_vals`: vector; containing the imposed dispacements;
- `disp_dofs`: DOFs interested by the Dirichlet BCs.

Returns a BoundaryConditions structure.
"""
function constructor_boundary_conditions(fixed_dofs, free_dofs, disp_fun::TF, disp_vals, disp_dofs, T=Float64) where TF

    return  BoundaryConditions{T, TF}(fixed_dofs, free_dofs, disp_fun, disp_vals, disp_dofs)
    
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
function constructor_configuration(mat::Material{TT, TJ}, geom::Geometry, ndofs, ext_forces::ExternalForces{TF}, bcs::BoundaryConditions{TT, TU}, T=Float64) where {TT, TF, TU, TJ}
    
    # displacement dofs
    disp_dofs = [i for i in 1:ndofs if mod1(i, 6)≤3]
    rot_dofs = [i for i in 1:ndofs if mod1(i, 6)>3]

    return Configuration{T, TJ, TF, TU}(mat, geom, ndofs, disp_dofs, rot_dofs, ext_forces, bcs)
    
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
function constructor_nodal_solution(Ktan, Ktan_free, ndofs, nfreedofs, T)

    return NodalSolution{T}(

    zeros(T, ndofs), #D
    zeros(T, ndofs), #Ḋ 
    zeros(T, ndofs), #D̈ 
    
    zeros(T, ndofs), #r
    Ktan, #Ktan_mat
    nonzeros(Ktan), #Ktan
    zeros(T, ndofs), #ΔD    
    zeros(T, ndofs), #temp vector

    zeros(T, nfreedofs), #r_free
    Ktan_free, #Ktan_free
    zeros(T, nfreedofs)) #ΔD_free

end

# Constructor of the structure containing the global matrices 
function constructor_global_matrices(I, J, ndofs, T=Float64)
    
    C_mat = sparse(I, J, zero(T))
    C = nonzeros(C_mat)
    
    M_mat = sparse(I, J, zero(T))
    M = nonzeros(M_mat)

    K_mat = sparse(I, J, zero(T))
    K =  nonzeros(K_mat)

    Tⁱⁿᵗ = zeros(T, ndofs)
    Tᵏ =  zeros(T, ndofs)
    Tᶜ = zeros(T, ndofs)
    Tᶜᵒⁿ = zeros(T, ndofs)

    return Matrices{T}(K_mat, K, C_mat, C, M_mat, M, Tⁱⁿᵗ, Tᵏ, Tᶜ, Tᶜᵒⁿ)
    
end 

# Constructor of the sparse matrices
function constructor_sparse_matrices!(beams, nodes, constraints, conf::Configuration{T}) where T

    ndofs =  conf.ndofs
    free_dofs = conf.bcs.free_dofs
    nfreedofs = length(free_dofs)

    I, J = sparsity(nodes, beams, constraints)

    Ktan = sparse(I, J, 0.)
    
    Ktan_free = Ktan[free_dofs, free_dofs]

    matrices = constructor_global_matrices(I, J, ndofs, T)
    nodes_sol = constructor_nodal_solution(Ktan, Ktan_free, ndofs, nfreedofs, T)

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
