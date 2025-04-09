# Compute contact gap, gradient, and Hessian for a spherical SDF.
@inline function contact_gap_beams(point, master_surface::SphereSurface, radius_beam)
    # Vector from the point to the center of the sphere
    aux = point - master_surface.center
    
    # Compute the norm (distance) from the point to the surface center
    norm_aux = norm(aux)
    
    # Gap is the distance to the surface minus the radius of the beam
    gₙ = norm_aux - master_surface.radius - radius_beam
    
    # Gradient of the gap with respect to the point's coordinates
    invnorm = 1 / norm_aux  # Inverse of the distance
    ∂gₙ∂x = invnorm * aux  # Gradient (direction of normal)
    
    # Hessian of the gap (second derivative)
    ∂²gₙ∂x² = invnorm * ID3 - (invnorm^3) * (aux * aux')
    
    return gₙ, ∂gₙ∂x, ∂²gₙ∂x²
end

# Determine if a point is in contact with a spherical SDF based on a threshold distance.
@inline function incontact_beams(point, master_surface::SphereSurface, radius_beam::Float64)

    # Define the threshold distance
    threshold_gap = radius_beam/4

    # Vector from the point to the center of the sphere
    aux = point - master_surface.center
    
    # Compute the distance to the surface and compare to the threshold
    gap = norm(aux) - master_surface.radius - radius_beam
    
    # If gₙ is less than or equal to the threshold gap, then the point is in contact
    return gap ≤ threshold_gap
end

# Quadratically regularize the penalty for contact calculations.
@inline function regularize_gap_penalty_beams(gap::Float64, radius_beam::Float64)
  
    # Define the threshold gap
    threshold_gap = radius_beam/4 

    # Define the regularization parameter as half of the threshold gap
    half_threshold_gap = threshold_gap / 2

    # If the gap is non-positive, use a simple quadratic penalty formulation
    if gap ≤ 0
        regularized_penalty = half_threshold_gap - gap  # Regularized penalty
        penalty_derivative = -one(Float64)  # Derivative of the penalty
        # Penalty energy: quadratic term and linear penalty for negative gap
        penalty_energy = gap^2 / 2 - half_threshold_gap * gap + (threshold_gap^2) / 6
    else
        # If gap is positive, compute the penalty for a smooth transition
        aux = (threshold_gap - half_threshold_gap) / (threshold_gap^2)
        regularized_penalty = aux * gap^2 - gap + half_threshold_gap  # Regularized penalty for positive gap
        penalty_derivative = 2 * aux * gap - 1  # Derivative of the penalty
        # Penalty energy for positive gaps
        penalty_energy = (threshold_gap - half_threshold_gap) / (3 * threshold_gap^2) * gap^3 - gap^2 / 2 + half_threshold_gap * gap - threshold_gap^2 / 6
    end
    
    return regularized_penalty, penalty_derivative, penalty_energy  # Return the penalty, its derivative, and penalty energy
end

# Smoothstep function for transition smoothing in contact calculations.
@inline function smoothstep_transition_beams(transition_value::Float64, position::Float64, radius_beam::Float64)

    upper_limit = radius_beam/4 

    # If the position is less than or equal to zero, return the transition value unchanged
    if position ≤ 0
        smoothed_value = transition_value
        smoothed_derivative = zero(Float64)
    else
        # Normalize the position within the range [0, upper_limit]
        normalized_position = (position - upper_limit) / upper_limit
        # Apply smoothstep function: cubic interpolation
        smoothed_value = normalized_position * normalized_position * (3 + 2 * normalized_position) * transition_value
        # Derivative of the smoothstep function
        smoothed_derivative = 6 / upper_limit * normalized_position * (1 + normalized_position) * transition_value
    end
    return smoothed_value, smoothed_derivative  # Return the smoothed value and its derivative
end