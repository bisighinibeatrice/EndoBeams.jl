@with_kw struct Params

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
    tcompt_max::Float64 = 1800
    stop_long_simulation::Bool = true 

    # HHT time stepping parameters
    α::Float64 = -0.05
    β::Float64 = 0.25*(1-α)^2
    γ::Float64 = 0.5*(1-2*α)

    # Solver tolerance and maximum number of iterations
    tol_res::Float64 = 1e-5
    tol_ΔD::Float64 = 1e-5
    max_it::Int = 10

    # Integration points
    nᴳ::Int = 3
    ωᴳ::Vec3{Float64} = Vec3(5/9, 8/9, 5/9)
    zᴳ::Vec3{Float64} = Vec3(-sqrt(3/5), 0, sqrt(3/5))

    #Constraint penalty
    kᶜᵒⁿ::Float64 = 1e3
    ηᶜᵒⁿ::Float64 = 1


end 

