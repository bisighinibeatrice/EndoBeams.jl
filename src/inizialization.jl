# Cleans the output folders from files of the precedent computation
function clean_folders(output_dir)

    if output_dir != ""
        dir = pwd()
        cd(output_dir)
        foreach(rm, filter(endswith(".vtp"), readdir()))
        foreach(rm, filter(endswith(".vtu"), readdir()))
        foreach(rm, filter(endswith(".vtk"), readdir()))
        foreach(rm, filter(endswith(".pvd"), readdir()))
        foreach(rm, filter(endswith(".txt"), readdir()))
        cd(dir)
    end   
end 

# Initializes simulation state, pre-allocates memory, and generates initial VTK files.
function setup_state_simulation(conf::BeamsConfiguration, inter::Union{Nothing, Interaction}, output_dir)
 
    # Allocate state variables for beams
    forcesⁿ = Forces(conf)
    forcesⁿ⁺¹ = deepcopy(forcesⁿ)
    matricesⁿ⁺¹, solutionⁿ⁺¹ = sparse_matrices_beams!(conf)
    energyⁿ⁺¹ = Energy()

    # Group beam state variables into a single structure
    state = SimulationState(forcesⁿ, forcesⁿ⁺¹, matricesⁿ⁺¹, solutionⁿ⁺¹, energyⁿ⁺¹)

    # Prepare VTK data for visualization
    vtkdata = VTKDataBeams(conf, output_dir)

    return state, VTKData(vtkdata)

end

# Initializes the simulation state for a beam-only configuration
function initialize_state_simulation!(conf::BeamsConfiguration, state::SimulationState, params::SimulationParams)
    
    # Assemble system matrices and initialize state variables
    assemble!(conf, state, params)
    if !isnothing(conf.constraints)
        assembly_constraints!(conf, state)
    end 

end
