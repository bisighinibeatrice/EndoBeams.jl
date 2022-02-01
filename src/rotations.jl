# Get skew symmetric matrix from angle
@inline function skew_skymmetric_matrix_from_vector(vec)

    return Mat33(0, vec[3], -vec[2], -vec[3], 0, vec[1], vec[2], -vec[1], 0)

end

# Get angle from skew symmetric matrix
@inline function get_angle_from_skew_skymmetric_matrix(SΘ)

    return Vec3(-SΘ[2,3], SΘ[1,3], -SΘ[1,2])

end

# Get the inverse of skew symmetric matrix from angle
function inverse_skew_skymmetric_matrix_from_angle(Θ::AbstractVector{T}) where T

    Θ_norm = norm(Θ)

    if Θ_norm < 10*eps(T)

        Tₛ⁻¹ = Mat33{T}(1,0,0,0,1,0,0,0,1)

    else
        SΘ = skew_skymmetric_matrix_from_vector(Θ)

        aux1 = Θ_norm/(2*tan(Θ_norm/2))

        aux2 = Θ/Θ_norm

        Tₛ⁻¹ = aux1*ID3 + (1-aux1)*(aux2*aux2') - SΘ/2 # in [Aguirre]: 15b
    end

    return Tₛ⁻¹

end

# Get Ts from angle
function get_Ts_from_angle(Θ::AbstractVector{T}) where T

    Θ_norm = norm(Θ)

    if Θ_norm < 10*eps(T)

        Tₛ⁻¹ = Mat33{T}(1,0,0,0,1,0,0,0,1)

    else
        SΘ = skew_skymmetric_matrix_from_vector(Θ)

        Tₛ⁻¹ = ID3 + (1-cos(Θ_norm))/(Θ_norm^2)*SΘ + (Θ_norm-sin(Θ_norm))/(Θ_norm^3)*(SΘ*SΘ) # in [Aguirre]: 15a
    end

    return Tₛ⁻¹

end

# Get angle from rotaion matrix (not computing the log)
function angle_from_rotation_matrix(R::AbstractMatrix{T}) where T

    if abs((tr(R)-1)/2) < 1
        norm_v = max(acos((tr(R)-1)/2), 10*eps(T))
    else
        norm_v = 10*eps(T)
    end

    n_v = (1/(2*sin(norm_v))) * Vec3(R[3,2]-R[2,3], R[1,3]-R[3,1], R[2,1]-R[1,2])

    return norm_v*n_v

end


# Rotate using Rodrigue's formula
function rotate_rod(a::AbstractVecOrMat{T}, Θ) where T
    
    Θ_norm = norm(Θ)
    if Θ_norm > 10*eps(T)
        K = skew_skymmetric_matrix_from_vector(Θ/Θ_norm)
        R = ID3 + sin(Θ_norm)*K + (1-cos(Θ_norm))*K^2
        return R*a
    else
        return a
    end

end



