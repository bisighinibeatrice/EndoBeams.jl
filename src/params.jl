@with_kw struct Params{T}

    # Options
    verbose::Bool = true
    record_timings::Bool = false
    output_dir::String = "output3D"

    # Time stepping
    ini_Δt::Float64 = 1e-2
    min_Δt::Float64 = 1e-10
    max_Δt::Float64 = 1e-2
    Δt_plot::Float64 = 1e-2
    tᵉⁿᵈ::Float64 = 1e-2
    accelerate_after_success_it::Int = 4
    stop_on_energy_threshold::Bool = false 
    energy_threshold::Float64 = 0.

    # HHT time stepping parameters
    α::T = -0.05
    β::T = 0.25*(1-α)^2
    γ::T = 0.5*(1-2*α)

    # Solver tolerance and maximum number of iterations
    tol_res::Float64 = 1e-5
    tol_ΔD::Float64 = 1e-5
    max_it::Int = 10

    # Integration points
    nᴳ::Int = 3
    ωᴳ::Vec3{T} = Vec3(5/9, 8/9, 5/9)
    zᴳ::Vec3{T} = Vec3(-sqrt(3/5), 0, sqrt(3/5))

    #Constraint penalty
    kᶜᵒⁿ::T = 1e3
    ηᶜᵒⁿ::T = 1


end 

