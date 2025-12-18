using Pkg; Pkg.activate(".")
using DelimitedFiles, ReadVTK, WriteVTK, LinearAlgebra, Test, Revise, EndoBeams, BenchmarkTools, StaticArrays, Dierckx, StructArrays, Parameters

# cd("stent")
include("utils/utils_IO.jl")
include("utils/utils_build_braided_stent_geometry.jl")
include("utils/utils_morping_cl.jl")
include("utils/utils_positioning.jl")

include("crimping.jl")
include("morph_centerline.jl")
include("positioning.jl")
include("deployment.jl") 

# --------------------------------------------
# Create stent
# --------------------------------------------

@unpack nbWires, rStent, rCrimpedStent, rWireSection, wireGap, lengthStent, nbTotalCells, braidingPattern = BraidedStent()
free_positions, connectivity = generate_braided_stent_geometry(nbWires, rStent, rWireSection, wireGap, lengthStent, nbTotalCells, braidingPattern)
constraints_connectivity = find_close_node_pairs_in_rings(free_positions)
free_positions = reduce(hcat, free_positions)'
connectivity = reduce(hcat, connectivity)'
constraints_connectivity = reduce(hcat, constraints_connectivity)'

# --------------------------------------------
# Crimping
# --------------------------------------------

output_dir_crimping = "stent_deployment/output3D/outputCrimping/"
if !isdir(output_dir_crimping) mkpath(output_dir_crimping) end

crimping(rStent, free_positions, connectivity, constraints_connectivity, output_dir_crimping)

# --------------------------------------------
# Geometrical morphing of the centerline
# --------------------------------------------

output_dir_morphing_cl = "stent_deployment/output3D/outputClMorphing/"
if !isdir(output_dir_morphing_cl) mkpath(output_dir_morphing_cl) end

filename_target_centerline = "stent_deployment/input/vessel_centerline.vtk"

# Morphing parameters
num_morphing_iterations = 100 # number of iterations for incremental morphing
deployment_offset_mm = 1.0 # offset along the target centerline in mm

deployment_origin_point = morph_stent_centerline_to_deployment!(free_positions, num_morphing_iterations, deployment_offset_mm, filename_target_centerline, output_dir_crimping, output_dir_morphing_cl)

# --------------------------------------------
# Positioning stent along centerline
# --------------------------------------------

output_dir_positioning = "stent_deployment/output3D/outputPositioning/"
if !isdir(output_dir_positioning) mkpath(output_dir_positioning) end

positioning(free_positions, connectivity, constraints_connectivity, num_morphing_iterations, deployment_origin_point, output_dir_crimping, output_dir_morphing_cl, output_dir_positioning)

# --------------------------------------------
# Deployment
# --------------------------------------------

output_dir_deployment = "stent_deployment/output3D/outputDeployment/"
if !isdir(output_dir_deployment) mkpath(output_dir_deployment) end

deployment(free_positions, connectivity, constraints_connectivity, deployment_origin_point, output_dir_crimping, output_dir_positioning, output_dir_deployment)

