@with_kw struct Params

    scale::Int = 10;  @assert scale == 2 || scale == 10
    ENERGY_STOP::Bool = false 
    SHOW_COMP_TIME::Bool = true
    SHOW_TIME_SECTIONS::Bool = false
    SAVE_NODES_VTK::Bool = true
    SAVE_ENERGY::Bool = false
    SAVE_INTERPOLATION_VTK::Bool = true
    SAVE_GP_VTK::Bool = false
    thisDirOutputPath::String = "output3D"

end 

@with_kw struct ParamsTest

    scale::Int = 2;  @assert scale == 2 || scale == 10
    ENERGY_STOP::Bool = false 
    SHOW_COMP_TIME::Bool = false
    SHOW_TIME_SECTIONS::Bool = false
    SAVE_ENERGY::Bool = false
    SAVE_NODES_VTK::Bool = false
    SAVE_INTERPOLATION_VTK::Bool = false
    SAVE_GP_VTK::Bool = false
    thisDirOutputPath::String = ""

end 
