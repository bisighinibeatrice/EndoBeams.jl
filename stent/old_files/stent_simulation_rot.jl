using DelimitedFiles
using Parameters
using LinearAlgebra
using Dierckx
using StaticArrays
# using SignedDistanceField
using GeometryBasics
using LinearAlgebra
using StructArrays
using GLMakie
using FileIO, MeshIO

# cd("stent")
include("utils_stent.jl")
include("crimping.jl")
include("positioning_cl.jl")
include("positioning.jl")
include("deployment.jl")

# -------------------------------------------------------------------------------------------
# Create stent
# -------------------------------------------------------------------------------------------

@unpack nbWires, rStent, rCrimpedStent, rWireSection, wireGap, lengthStent, nbTotalCells, braidingPattern = BraidedStent()
initial_positions_stent, connectivity_stent = compute_bs_geom_given_nbTotalCells(nbWires, rStent, rWireSection, wireGap, lengthStent, nbTotalCells, braidingPattern)
initial_positions_stent_mat = reshape(reinterpret(Float64, initial_positions_stent), (3, length(initial_positions_stent)))'

output_dir_crimping = "stent/output3D/outputCrimping3D/"
crimping(rStent, rCrimpedStent, initial_positions_stent, connectivity_stent, output_dir_crimping)

output_dir_positioning_cl = "stent/output3D/outputPositioningCl3D/"
output_dir_positioning = "stent/output3D/outputPositioning3D/"
output_dir_deployment = "stent/output3D/outputDeployment3D/"

filename_cl = "stent/input/cl_model.vtk"
filename_surf = "stent/input/model.stl"
filename_sdf = "stent/input/sdf.vtk"

nb_iterations = 20
deploy_pos =  1

# -------------------------------------------------------------------------------------------
# Reorientation
# -------------------------------------------------------------------------------------------

# cl_T = read_vtk_centerline(filename_cl)
# v_η = vcat(0, cumsum(norm.(diff(cl_T))))
# spline = Vec3(Dierckx.Spline1D(v_η, getindex.(cl_T, i)) for i in 1:3)

# crimped_positions_stent = initial_positions_stent + read_ics_vec(readdlm(output_dir_crimping * "u.txt"))
# pos_cl_0 = get_centerline_stent(crimped_positions_stent)

# # get final centerline
# cl_η = vcat(0, cumsum(norm.(diff(pos_cl_0))))
# npointsCl = length(cl_η)
# pos_cl_T = zeros(Vec3, npointsCl)
# for i in 1:npointsCl
#     pos_cl_T[i] = Dierckx.evaluate.(spline, deploy_pos + cl_η[i])
# end 
# set_origin!(crimped_positions_stent, pos_cl_0[1]-pos_cl_T[1])   
# pos_cl_0 =  get_centerline_stent(crimped_positions_stent, pos_cl_T[1])

# surface = load(filename_surf)
# positions = decompose(Point3{Float32}, surface)
# trisconn = decompose(TriangleFace{Int}, surface)
# positions_rot = rotate(positions, pos_cl_T[1], pos_cl_0[end]-pos_cl_0[1], pos_cl_T[end]-pos_cl_T[1])
# new_mesh =   MeshIO.Mesh(positions_rot, trisconn)

filename_surf = "stent/input/rot_model.stl"
# save(filename_surf, new_mesh)

# cl__rot = rotate(cl_T, pos_cl_T[1], pos_cl_0[end]-pos_cl_0[1], pos_cl_T[end]-pos_cl_T[1])

filename_cl = "stent/input/cl_rot_model.vtk"
# write_vtk_configuration(filename_cl, cl__rot, [])

# -------------------------------------------------------------------------------------------
# Geometrical positioning centerline
# -------------------------------------------------------------------------------------------

initial_cl_stent = positioning_cl(initial_positions_stent, connectivity_stent, nb_iterations, deploy_pos, filename_cl, output_dir_crimping, output_dir_positioning_cl)

# -------------------------------------------------------------------------------------------
# Physical positioning stent
# -------------------------------------------------------------------------------------------

positioning(initial_positions_stent, connectivity_stent, nb_iterations, initial_cl_stent[1], output_dir_crimping, output_dir_positioning_cl, output_dir_positioning)

crimped_positions_stent = initial_positions_stent .+  read_ics_vec(readdlm(output_dir_crimping * "u.txt"))
positions_cl =  get_centerline_stent(crimped_positions_stent)
set_origin!(crimped_positions_stent, positions_cl[1]-initial_cl_stent[1])
open(output_dir_positioning * "/crimped.txt", "w") do io
    writedlm(io, crimped_positions_stent)
end

disp = read_ics_vec(readdlm(output_dir_positioning * "u.txt"))
write_vtk_configuration("stent/output3D/positioned.vtk", crimped_positions_stent + disp, connectivity_stent)

# -------------------------------------------------------------------------------------------
# Deployment
# -------------------------------------------------------------------------------------------

# if isfile(filename_sdf)
#     rm(filename_sdf)
# end
# s, x, y, z = sdfgen(filename_surf, 0.1; padding = 1, acceleration=:KDTree)
# write_vtk_general_structured_mesh(filename_sdf, -s, 0.1, x, y, z)

deployment(initial_positions_stent_mat, initial_positions_stent, connectivity_stent, filename_sdf, output_dir_crimping, output_dir_positioning, output_dir_deployment)

disp = read_ics_vec(readdlm(output_dir_deployment * "u.txt"))
write_vtk_configuration("stent/output3D/deployed.vtk", initial_positions_stent + disp, connectivity_stent)

