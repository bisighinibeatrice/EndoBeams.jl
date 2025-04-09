#----------------------------------
# STRUCTURES FOR VISUALISATION
#----------------------------------

# Structure to store VTK data
# The structure holds data for VTK visualization if available
struct VTKData
   
    # VTK data for beams, can be `Nothing` if beams are not part of the simulation
    vtkdata_beams::Union{Nothing, VTKDataBeams} 
    
end

#----------------------------------
# FUNCTION TO WRITE VTK DATA
#----------------------------------

function write_VTK(write_counter::Int, step::Int, t::Float64, conf::BeamsConfiguration, vtkdata::VTKData)

    write_VTK_beams(write_counter, step, t, conf, vtkdata.vtkdata_beams)
    
end
