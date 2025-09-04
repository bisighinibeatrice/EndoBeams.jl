using Pkg; Pkg.activate(".")
using DelimitedFiles, ReadVTK, WriteVTK, LinearAlgebra, Test, Revise, EndoBeams, BenchmarkTools, StaticArrays, Dierckx, StructArrays, Parameters

# cd("stent")
include("utils/utils_IO.jl")
include("utils/utils_stent.jl")
include("utils/utils_build_braided_stent_geometry.jl")
include("utils/utils_morping_cl.jl")
include("utils/utils_positioning.jl")

include("crimping.jl")
include("morph_centerline.jl")
include("positioning.jl")

# -------------------------------------------------------------------------------------------
# Create stent
# -------------------------------------------------------------------------------------------

rStent = 2.3
positions_stent = readdlm("positions_stent.txt")
connectivity_stent = readdlm("connectivity_stent.txt", Int)
constraints_connectivity_stent = readdlm("constraints_stent.txt", Int)

# -------------------------------------------------------------------------------------------
# Crimping
# -------------------------------------------------------------------------------------------

output_dir_crimping = "stent/output3D/outputCrimping3D/"
if !isdir(output_dir_crimping) mkpath(output_dir_crimping) end

crimping(rStent, positions_stent, connectivity_stent, constraints_connectivity_stent, output_dir_crimping)

# -------------------------------------------------------------------------------------------
# Geometrical morphing of the centerline
# -------------------------------------------------------------------------------------------

output_dir_positioning_cl = "stent/output3D/outputClMorphing3D/"
if !isdir(output_dir_positioning_cl) mkpath(output_dir_positioning_cl) end

filename_cl = "stent/input/cl_model.vtk"
nb_iterations = 100
deploy_pos =  1

initial_cl_stent = morph_stent_centerline_to_deployment!(positions_stent, connectivity_stent, nb_iterations, deploy_pos, filename_cl, output_dir_crimping, output_dir_positioning_cl)

# -------------------------------------------------------------------------------------------
# Positioning stent along centerline
# -------------------------------------------------------------------------------------------

output_dir_positioning = "stent/output3D/outputPositioning3D/"
if !isdir(output_dir_positioning) mkpath(output_dir_positioning) end

positioning(positions_stent, connectivity_stent, constraints_connectivity_stent, nb_iterations, initial_cl_stent[1], output_dir_crimping, output_dir_positioning_cl, output_dir_positioning)
