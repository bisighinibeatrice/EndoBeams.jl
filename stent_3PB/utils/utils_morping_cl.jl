"""
Compute segment vectors and their normalized directions from a list of positions.

# Arguments
- `positions::Vector{Vec3}`: Sequence of 3D points along a path.
- `ref_vec::Vec3`: Reference direction vector for the first segment (default: (0, 0, 1)).

# Returns
- `segments::Vector{Vec3}`: Vectors between successive positions.
- `normalized_segments::Vector{Vec3}`: Unit vectors of the above segments.
"""
function compute_segments_and_directions(positions::Vector{Vec3}, ref_vec::Vec3 = Vec3(0, 0, 1))
    n = length(positions)

    segments = Vector{Vec3}(undef, n)
    normalized_segments = Vector{Vec3}(undef, n)

    segments[1] = ref_vec

    for i in 2:n
        segments[i] = positions[i] - positions[i - 1]
    end

    for i in 1:n
        normalized_segments[i] = segments[i] / norm(segments[i])
    end

    return segments, normalized_segments
end

"""
Compute the rotation axis between two unit vectors `a` and `b`.

# Returns
- A unit vector representing the axis of rotation from `a` to `b`.
"""
function compute_rotation_axis(a::Vec3, b::Vec3)
    @assert isapprox(norm(a), 1.0; atol=1e-8) "Vector `a` must be normalized."
    @assert isapprox(norm(b), 1.0; atol=1e-8) "Vector `b` must be normalized."

    axis = cross(a, b)
    if iszero(norm(axis))
        return Vec3(0.0, 0.0, 1.0)  # Arbitrary fallback axis
    else
        return axis / norm(axis)
    end
end

"""
Compute the rotation angle (in radians) between two unit vectors.

# Returns
- Rotation angle in radians.
"""
function compute_rotation_angle(a::Vec3, b::Vec3)
    @assert isapprox(norm(a), 1.0; atol=1e-8) "Vector `a` must be normalized."
    @assert isapprox(norm(b), 1.0; atol=1e-8) "Vector `b` must be normalized."

    return acos(clamp(dot(a, b), -1.0, 1.0))  # Clamp avoids domain errors
end

"""
Compute a rotation matrix using Rodrigues' rotation formula.

# Arguments
- `θ::Float64`: Rotation angle (in radians).
- `axis::Vec3`: Unit vector axis of rotation.

# Returns
- `Mat33`: Rotation matrix.
"""
function rodrigues_rotation_matrix(θ::Float64, axis::Vec3)
    c = cos(θ)
    s = sin(θ)
    t = 1 - c

    ux, uy, uz = axis
    uxx, uyy, uzz = ux^2, uy^2, uz^2
    uxy, uxz, uyz = ux*uy, ux*uz, uy*uz

    return Mat33(
        c + uxx * t,    uxy * t + uz * s,  uxz * t - uy * s,
        uxy * t - uz * s, c + uyy * t,    uyz * t + ux * s,
        uxz * t + uy * s, uyz * t - ux * s, c + uzz * t
    )
end

"""
Rotate point `p` around center `c` using rotation matrix `M`.

# Arguments
- `p::Vec3`: Point to rotate.
- `c::Vec3`: Rotation center.
- `M::Mat33`: Rotation matrix.

# Returns
- `Vec3`: Rotated point.
"""
function rotate_point_around_center(p::Vec3, c::Vec3, M::Mat33)
    return M * (p - c) + c
end
