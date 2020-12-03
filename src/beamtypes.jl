# The structs defining the beam types should be put there

struct Node{T}
	coords::Vec3{T}
	# and so on
end


struct Beam{T}
	n1idx::Int
	n2idx::Int
	initial_length::T
	# and so on
end