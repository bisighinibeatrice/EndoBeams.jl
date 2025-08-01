using Pkg; Pkg.activate(".")
using DelimitedFiles, ReadVTK, WriteVTK, LinearAlgebra, Test, Revise, EndoBeams, BenchmarkTools, StaticArrays, Dierckx, StructArrays, Parameters

# cd("stent")
include("utils_stent.jl")
include("crimping.jl")

# -------------------------------------------------------------------------------------------
# Create stent
# -------------------------------------------------------------------------------------------

rStent = 2.3
positions_stent = readdlm("stent/input/positions_stent.txt")
connectivity_stent = readdlm("stent/input/connectivity_stent.txt", Int)
constraints_connectivity_stent = readdlm("stent/input/constraints_stent.txt", Int)

# -------------------------------------------------------------------------------------------
# Crimping
# -------------------------------------------------------------------------------------------

output_dir_crimping = "stent/output3D/outputCrimping3D/"
if !isdir(output_dir_crimping) mkpath(output_dir_crimping) end
crimping(rStent, positions_stent, connectivity_stent, constraints_connectivity_stent, output_dir_crimping)