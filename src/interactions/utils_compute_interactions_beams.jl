# Compute contact gap, gradient, and Hessian for a spherical SDF.
@inline function contact_gap_beams(point, master_surface::SphereSurface, radius_beam::Float64)
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
    ḡₙ = radius_beam/4

    # Vector from the point to the center of the sphere
    aux = point - master_surface.center
    
    # Compute the distance to the surface and compare to the threshold
    gₙ = norm(aux) - master_surface.radius - radius_beam
    
    # If gₙ is less than or equal to the threshold gap, then the point is in contact
    return gₙ ≤ ḡₙ
end

# Compute contact gap, gradient, and Hessian for a planar SDF.
@inline function contact_gap_beams(point, master_surface::PlaneSurface, radius_beam::Float64)
    # Extract plane data
    x₀ = master_surface.point
    n  = master_surface.normal   # must be unit
    
    # Gap: signed distance to plane minus beam radius
    gₙ = dot(n, point - x₀) - radius_beam

    # Gradient is constant for a plane
    ∂gₙ∂x = n

    # Hessian is identically zero
    ∂²gₙ∂x² = zeros(3,3)

    return gₙ, ∂gₙ∂x, ∂²gₙ∂x²
end

# Determine if point is in contact with planar SDF
@inline function incontact_beams(point, master_surface::PlaneSurface, radius_beam::Float64)
    # Threshold (same style as sphere version)
    ḡₙ = radius_beam/4

    x₀ = master_surface.point
    n  = master_surface.normal

    # Gap
    gₙ = dot(n, point - x₀) - radius_beam

    # Contact when gap <= threshold
    return gₙ ≤ ḡₙ
end

# Check if a point is within the specified domain.
@inline function isinside(point, dom)
    l_x = point[1] - dom[1]  
    flag_lx = l_x >= 0 && l_x <= (dom[2] - dom[1]) 
    l_y = point[2] - dom[3]
    flag_ly = l_y >= 0 && l_y <= (dom[4] - dom[3])
    l_z = point[3] - dom[5] 
    flag_lz = l_z >= 0 && l_z <= (dom[6] - dom[5])
    return flag_lx && flag_ly && flag_lz
end

# Compute contact gap, gradient, and Hessian for a discrete SDF.
@inline function contact_gap_beams(point, sdf::DiscreteSignedDistanceField, radius_beam::Float64)
    
    sitp = sdf.sitp
    gₙ = sitp(point...)
    ∂gₙ∂x = Interpolations.gradient(sitp, point...)
    ∂²gₙ∂x² = Interpolations.hessian(sitp, point...)
    
    # Normalize gradient for continuous contact calculations.
    nn = dot(∂gₙ∂x, ∂gₙ∂x)
    nmaginv = 1 / sqrt(nn)
    ∂gₙ∂x = ∂gₙ∂x * nmaginv 
    ∂²gₙ∂x² = ∂²gₙ∂x² * nmaginv * (ID3 - (∂gₙ∂x * ∂gₙ∂x') / nn)
    
    return gₙ - radius_beam, ∂gₙ∂x, ∂²gₙ∂x²
end 

# Check if a point is in contact with a discrete SDF.
@inline function incontact_beams(point, sdf::DiscreteSignedDistanceField, radius_beam::Float64)

    # Define the threshold distance
    ḡₙ = radius_beam/4

    # Check if the point is inside the SDF domain and compute the gap
    return isinside(point, sdf.dom) && sdf.sitp(point...) - radius_beam ≤ ḡₙ
    
end 
