using EndoBeams
using Test

const T = Float32

@testset "EndoBeams.jl" begin

    println("Start 1st test")
    include("test_angle.jl")
    println("End 1st test")

    # println("Start 2nd test")
    # include("test_beam_contact_sphere.jl")
    # println("End 2nd test")

    # println("Start 3rd test")
	# include("test1_ring_contact_plane.jl")
    # println("End 3rd test")

    # println("Start 4th test")
    # include("test2_ring_contact_plane.jl")
    # println("End 4th test")

    # println("Start 5th test")
    # include("test_ring_crimping.jl")
    # println("End 5th test")

end

