# Computes the skew-symmetric matrix of a 3D vector.
@inline function skew(vec)
    
    # Create the 3x3 skew-symmetric matrix
    return @SMatrix [  0       -vec[3]  vec[2] ;    # First row:  [0, -z,  y]
                      vec[3]   0      -vec[1];    # Second row: [z,  0,  -x]
                      -vec[2]  vec[1]   0 ];    # Third row:  [-y, x,  0]
end

# Computes the rotation matrix using Rodrigues' rotation formula.
@inline function rotation_matrix(Θ::AbstractVecOrMat{Float64})

    Θ_norm = norm(Θ)  # Magnitude of the rotation vector, which is also the angle of rotation (θ)
    
    # If the rotation angle is significantly small (close to 0), return the identity matrix
    if Θ_norm > 10 * eps(Float64)
        sinΘ = sin(Θ_norm)  # Sine of the rotation angle
        cosΘ = cos(Θ_norm)  # Cosine of the rotation angle
        
        # Skew-symmetric matrix of Θ, which represents the cross product operation
        SΘ = skew(Θ)
        
        # Rodrigues' rotation formula for computing the rotation matrix
        R = ID3 + sinΘ / Θ_norm * SΘ + (1 - cosΘ) / Θ_norm^2 * SΘ * SΘ
        
        return R
    else
        # For very small angles (where the rotation is close to zero), return the identity matrix
        return SMatrix{3,3,Float64,9}(1, 0, 0, 0, 1, 0, 0, 0, 1)
    end

end

# Computes the local rotation matrix and other auxiliary matrices and parameters.
@inline function local_Rₑ_and_aux(x₁, x₂, R₁, R₂, Rₑ⁰E₂, lₙ)
    # Compute the normalized direction vector from x₁ to x₂
    v₁ = (x₂ - x₁) / lₙ

    # Compute the projections of R₁ and R₂ onto Rₑ⁰E₂
    p₁ = R₁ * Rₑ⁰E₂
    p₂ = R₂ * Rₑ⁰E₂

    # Average the projections
    p = (p₁ + p₂) / 2

    # Compute the cross product of v₁ and p to get the third basis vector v₃
    v₃ = cross(v₁, p)
    v₃ = v₃ / norm(v₃)

    # Compute the second basis vector v₂ using the cross product
    v₂ = cross(v₃, v₁)

    # Assemble the rotation matrix Rₑ
    Rₑ = [v₁ v₂ v₃]

    # Define auxiliary matrices
    ru₁ = -v₁'  # Ru1 matrix
    ru₂ = v₁'   # Ru2 matrix

    # Extract 2x2 part of Rₑ matrix
    Rₑ¹²ᵀ = Rₑ[:, SOneTo(2)]'

    # Compute various auxiliary quantities
    q₁, q₂ = Rₑ¹²ᵀ * p
    q¹¹, q¹² = Rₑ¹²ᵀ * p₁
    q²¹, q²² = Rₑ¹²ᵀ * p₂

    # Compute eta and other parameters
    η = q₁ / q₂
    η¹¹ = q¹¹ / q₂
    η¹² = q¹² / q₂
    η²¹ = q²¹ / q₂
    η²² = q²² / q₂

    # Compute the matrix Gᵀu₁ and Gᵀu₂ for the rotations
    Gᵀu₁ = @SMatrix [0 0 η / lₙ; 0 0 1 / lₙ; 0 -1 / lₙ 0]
    Gᵀu₂ = -Gᵀu₁

    # Compute GᵀΘ₁ and GᵀΘ₂ matrices for the rotations
    GᵀΘ₁ = @SMatrix [η¹² / 2 -η¹¹ / 2 0; 0 0 0; 0 0 0]
    GᵀΘ₂ = @SMatrix [η²² / 2 -η²¹ / 2 0; 0 0 0; 0 0 0]

    # Compute the matrix D₃
    D₃ = (ID3 - v₁ * v₁') / lₙ

    return Rₑ, ru₁, ru₂, η, Gᵀu₁, GᵀΘ₁, Gᵀu₂, GᵀΘ₂, D₃
end

# Converts a rotation matrix to the corresponding angle (and axis) using Rodrigues' formula.
@inline function toangle(R::AbstractMatrix{Float64})

    # Check if the rotation matrix is close to the identity matrix (no rotation)
    if abs((tr(R) - 1) / 2) < 1
        norm_v = max(acos((tr(R) - 1) / 2), 10 * eps(Float64))
    else
        norm_v = 10 * eps(Float64)  # Small angle case, return a very small value
    end

    # Calculate the unit vector representing the axis of rotation
    n_v = (1 / (2 * sin(norm_v))) * @SVector [R[3, 2] - R[2, 3], R[1, 3] - R[3, 1], R[2, 1] - R[1, 2]]

    # Return the axis scaled by the rotation angle
    return norm_v * n_v
end

# Computes the inverse of the skew-symmetric transformation matrix.
@inline function Tₛ⁻¹(Θ::AbstractVector{Float64})
    Θ_norm = norm(Θ)  # Compute the magnitude of the rotation vector

    if Θ_norm < 10 * eps(Float64)
        # For small Θ, use the identity matrix as an approximation
        return SMatrix{3,3,Float64,9}(1, 0, 0, 0, 1, 0, 0, 0, 1)
    else
        SΘ = skew(Θ)  # Compute the skew-symmetric matrix of Θ
        sinΘ2, cosΘ2 = sincos(Θ_norm / 2)  # Compute sin and cos of half the angle

        # Compute the inverse transformation matrix using the derived formula
        Tₛ⁻¹ = ID3 - SΘ / 2 + ((sinΘ2 - Θ_norm / 2 * cosΘ2) / (sinΘ2 * Θ_norm^2)) * SΘ * SΘ

        return Tₛ⁻¹
    end
end

# Computes the transformation matrix Tₛ from a given rotation vector Θ.
@inline function Tₛ(Θ::AbstractVector{Float64})
    Θ_norm = norm(Θ)  # Compute the magnitude of the rotation vector

    if Θ_norm < 10 * eps(Float64)
        # For very small Θ, return the identity matrix as an approximation
        return SMatrix{3,3,Float64,9}(1, 0, 0, 0, 1, 0, 0, 0, 1)
    else
        SΘ = skew(Θ)  # Compute the skew-symmetric matrix of Θ

        # Compute the transformation matrix using the derived formula
        Tₛ = ID3 + (2 * (sin(Θ_norm / 2) / Θ_norm)^2) * SΘ +
                  ((1 - sin(Θ_norm) / Θ_norm) / Θ_norm^2) * (SΘ * SΘ)

        return Tₛ
    end
end

# Computes the Kᵥ matrix used in rotational kinematics.
@inline function compute_Kᵥ(Θ::AbstractVector{Float64}, v)
    Θnorm = norm(Θ)  # Compute the magnitude of the rotation vector

    if Θnorm < 10 * eps(Float64)
        # If rotation angle is very small, approximate Kᵥ using the skew-symmetric matrix
        return skew(v) / 2
    else
        # Compute trigonometric terms
        sinΘ, cosΘ = sincos(Θnorm)
        sinΘ2 = sin(Θnorm / 2)

        # Compute auxiliary term used in scaling factors
        aux = (2 * sinΘ2 / Θnorm)^2

        # Compute unit rotation axis (direction of Θ)
        u = Θ / Θnorm

        # Compute vector projections and cross-products
        uvᶜ = cross(u, v)   # Cross product of u and v
        uvᵈ = dot(u, v)     # Dot product of u and v

        # Compute outer products for matrix operations
        UU = u * u'  # Outer product of u with itself
        VU = v * u'  # Outer product of v with u
        UV = u * v'  # Outer product of u with v

        # Compute the transformation matrix Kᵥ using derived formula
        Kᵥ =  (cosΘ - sinΘ / Θnorm) / Θnorm * (VU - uvᵈ * UU) +
              (1 - sinΘ / Θnorm) / Θnorm * (UV - 2 * uvᵈ * UU + uvᵈ * ID3) -
              (sinΘ / Θnorm - aux) * (uvᶜ * u') +
              aux * skew(v) / 2

        return Kᵥ
    end
end

# Computes the η and μ coefficients used in rotational kinematics.
@inline function compute_η_μ(Θ̄::AbstractVector{Float64})
    Θ = norm(Θ̄)  # Compute the magnitude of the rotation vector

    if Θ < 10 * eps(Float64)
        # Small-angle approximation using Taylor series expansions
        η = Float64(1 / 12)
        μ = Float64(1 / 360)
    else
        # Compute trigonometric terms
        sinΘ, cosΘ = sincos(Θ)
        sinΘ2 = sin(Θ / 2)

        # Compute η and μ using derived mathematical expressions
        η = ((2 * sinΘ) - Θ * (1 + cosΘ)) / (2 * Θ^2 * sinΘ)
        μ = (Θ * (Θ + sinΘ) - 8 * sinΘ2^2) / (4 * Θ^4 * sinΘ2^2)
    end

    return η, μ
end

# Computes the correction matrix K̄ₕ used in rotational kinematics.
@inline function compute_K̄ₕ(Θ̅, M̄, Tₛ⁻¹Θ̅, η, μ)
    # Compute auxiliary matrices
    Θ̅M̄ᵀ = Θ̅ * M̄'    # Outer product of Θ̅ and M̄
    M̄Θ̅ᵀ = Θ̅M̄ᵀ'      # Transpose to reverse multiplication order

    SΘ̅ = skew(Θ̅)      # Skew-symmetric matrix of Θ̅
    SM̄ = skew(M̄)      # Skew-symmetric matrix of M̄

    # Compute K̄ₕ using the given mathematical formulation
    K̄ₕ = ( η * (Θ̅M̄ᵀ - 2 * M̄Θ̅ᵀ + dot(Θ̅, M̄) * ID3) + 
            μ * (SΘ̅ * SΘ̅ * M̄Θ̅ᵀ) - 
            SM̄ / 2 ) * Tₛ⁻¹Θ̅

    return K̄ₕ
end

# Computes the P matrices used in structural or kinematic computations.
@inline function Pmatrices(N₁, N₂, N₃, N₄, N₅, N₆, lₙ, η, η₁₁, η₁₂, η₂₁, η₂₂)

    # First set of P matrices
    P₁P¹ = @SMatrix [ 0  0  0; 
                      0  (N₃+N₄)/lₙ  0; 
                      0  0  (N₃+N₄)/lₙ ]

    P₁P² = @SMatrix [ 0  0  0; 
                      0  0  N₃; 
                      0  -N₃  0 ]

    P₁P³ = -P₁P¹  # Symmetric relationship
    P₁P⁴ = @SMatrix [ 0  0  0; 
                      0  0  N₄; 
                      0  -N₄  0 ]

    # Second set of P matrices
    P₂P¹ = @SMatrix [ 0  0  -η*(N₁+N₂)/lₙ; 
                      0  0  -(N₅+N₆)/lₙ; 
                      0  (N₅+N₆)/lₙ  0 ]

    P₂P² = @SMatrix [ -(N₂*η₁₂)/2 - N₁*(η₁₂/2 - 1)  (η₁₁*(N₁+N₂))/2  0; 
                       0  N₅  0; 
                       0  0  N₆ ]

    P₂P³ = -P₂P¹  # Symmetric relationship
    P₂P⁴ = @SMatrix [ -(N₁*η₂₂)/2 - N₂*(η₂₂/2 - 1)  (η₂₁*(N₁+N₂))/2  0; 
                       0  N₅  0; 
                       0  0  N₆ ]

    return P₁P¹, P₁P², P₁P³, P₁P⁴, P₂P¹, P₂P², P₂P³, P₂P⁴
end

