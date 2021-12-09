# Get skew symmetric matrix from angle
@inline function get_skew_skymmetric_matrix_from_vector(vec)

    return Mat33(0, vec[3], -vec[2], -vec[3], 0, vec[1], vec[2], -vec[1], 0)

end

# Get angle from skew symmetric matrix
@inline function get_angle_from_skew_skymmetric_matrix(Stheta)

    return Vec3(-Stheta[2,3], Stheta[1,3], -Stheta[1,2])

end

# Get the inverse of skew symmetric matrix from angle
function get_inverse_skew_skymmetric_matrix_from_angle(theta::AbstractVector{T}) where T

    theta_norm = norm(theta)

    if theta_norm < 10*eps(T)

        Tsinv = ID3

    else
        Stheta = get_skew_skymmetric_matrix_from_vector(theta)

        aux1 = theta_norm/(2*tan(theta_norm/2))

        aux2 = theta/theta_norm

        Tsinv = aux1*ID3 + (1-aux1)*(aux2*aux2') - Stheta/2 # in [Aguirre]: 15b
    end

    return Tsinv

end

# Get Ts from angle
function get_Ts_from_angle(theta::AbstractVector{T}) where T

    theta_norm = norm(theta)

    if theta_norm < 10*eps(T)

        Tsinv = ID3

    else
        Stheta = get_skew_skymmetric_matrix_from_vector(theta)

        Tsinv = ID3 + (1-cos(theta_norm))/(theta_norm^2)*Stheta + (theta_norm-sin(theta_norm))/(theta_norm^3)*(Stheta*Stheta) # in [Aguirre]: 15a
    end

    return Tsinv

end

# Get angle from rotaion matrix (not computing the log)
function get_angle_from_rotation_matrix(R::AbstractMatrix{T}) where T

    if abs((tr(R)-1)/2) < 1
        norm_v = max(acos((tr(R)-1)/2), 10*eps(T))
    else
        norm_v = 10*eps(T)
    end

    n_v = (1/(2*sin(norm_v))) * Vec3(R[3,2]-R[2,3], R[1,3]-R[3,1], R[2,1]-R[1,2])

    return norm_v*n_v

end


# Rotate using Rodrigue's formula
function rotate_rod(a::AbstractVecOrMat{T}, theta) where T
    
    theta_norm = norm(theta)
    if theta_norm > 10*eps(T)
        K = get_skew_skymmetric_matrix_from_vector(theta/theta_norm)
        R = ID3 + sin(theta_norm)*K + (1-cos(theta_norm))*K^2
        return R*a
    else
        return a
    end

end