@with_kw struct SimulationParams

    # Options
    verbose::Bool = true
    record_timings::Bool = false
    output_dir::String = "output3D"

    # Time stepping
    initial_timestep::Float64 = 1e-2
    min_timestep::Float64 = 1e-10
    max_timestep::Float64 = 1e-2
    output_timestep::Float64 = 1e-4
    simulation_end_time::Float64 = 1e-2
    accelerate_after_success_it::Int = 4
    stop_on_energy_threshold::Bool = false 
    energy_threshold::Float64 = 0.
    tcompt_max::Float64 = 1800
    stop_long_simulation::Bool = false 

    # HHT time stepping parameters
    α::Float64 
    β::Float64
    γ::Float64

    # Solver tolerance and maximum number of iterations
    tolerance_residual::Float64 = 1e-5
    tolerance_displacement::Float64 = 1e-5
    max_iterations::Int = 20
    
    # Integration points for beams
    nᴳ_beams::Int = 3
    ωᴳ_beams::Vec3{Float64} = Vec3(5/9, 8/9, 5/9)
    zᴳ_beams::Vec3{Float64} = Vec3(-sqrt(3/5), 0, sqrt(3/5))

end 

