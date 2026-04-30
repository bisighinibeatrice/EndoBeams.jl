using Pkg; Pkg.activate(".")
using DelimitedFiles, ReadVTK, WriteVTK, LinearAlgebra, Test, Revise, EndoBeams, BenchmarkTools, StaticArrays, Dierckx, StructArrays, Parameters

# cd("stent")
include("utils/utils_IO.jl")
include("utils/utils_build_braided_stent_geometry.jl")
include("utils/utils_morping_cl.jl")
include("utils/utils_positioning.jl")

include("crimping.jl")

# --------------------------------------------
# Create stent
# --------------------------------------------

free_positions = readdlm("stent_deployment/input/stent_positions.txt", ',')
connectivity = Int.(readdlm("stent_deployment/input/stent_connectivity.txt", ','))
# free_positions = reduce(hcat, free_positions)'
# connectivity = reduce(hcat, connectivity)'
rStent = 2

# --------------------------------------------
# Crimping
# --------------------------------------------

output_dir_crimping = "stent_deployment/output3D/outputCrimping/"
if !isdir(output_dir_crimping) mkpath(output_dir_crimping) end

crimping_ring(rStent, free_positions, connectivity, output_dir_crimping)