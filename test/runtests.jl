using EndoBeams
using Test

@testset "EndoBeams.jl" begin

    include("test_angle.jl")
    include("test_beam_contact_sphere.jl")
	include("test1_ring_contact_plane.jl")
    include("test2_ring_contact_plane.jl")
    include("test_ring_crimping.jl")

end

