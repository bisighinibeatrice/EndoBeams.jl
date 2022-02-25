
#----------------------------------
# STRUCTURES
#----------------------------------

"""
    sdf = struct SDF_Plane_z{T}

Constructor of the structure containing the properties of the analytical SDF of a z-normal plane:
- `r`: beam radius;
- `z0`: plane position along z.

Returns a structure containing the information of the created sdf. 
"""
struct SDF_Plane_z{T}
    
    r::T  # beams radius
    z0::T # plane position 
    
end 

"""
    sdf = struct SDF_Plane_y{T}

Constructor of the structure containing the properties  of the analytical SDF of a y-normal plane.
- `r`: beam radius;
- `y0`: plane position along y.

Returns a structure containing the information of the created sdf.  
"""
struct SDF_Plane_y{T}
    
    r::T  # beams radius
    y0::T # plane position 
    
end 

"""
    sdf = struct SDF_Sphere{T}

Constructor of the structure containing the properties  of the analytical SDF of a sphere.
- `r`: beam radius;
- `R`: sphere radius;
- `x0`: sphere centre position along x;
- `y0`: sphere centre position along y;
- `z0`: sphere centre position along z.
    
Returns a structure containing the information of the created sdf. 
"""
struct SDF_Sphere{T}
    
    r::T # beams radius
    R::T # sphere radius
    x0::T # sphere position
    y0::T # sphere position
    z0::T # sphere position
    
end 

"""
    sdf = struct SDF_Cylinder{T}

Constructor of the structure containing the properties of the analytical SDF of an infinite cylinder oriented along z.
- `r`: beam radius;
- `R`: cylinder radius.
    
Returns a structure containing the information of the created sdf. 
"""
struct SDF_Cylinder{T}
    
    r::T # beams radius
    R::T # cylinder radius
    
    
end  

# Properties of a discrete SDF
struct SDF_Discrete{T, F}
    
    r::T
    sitp::F
    dom::NTuple{6, Float64}
    dx::T
    dy::T
    dz::T
    
end 

#----------------------------------
# CONSTRUCTOR DISCRETE SDF
#----------------------------------
"""
    sdf = constructor_discrete_sdf(filename, radius, inside,  T=Float64)

Constructor of the discrete SDF from a vtk file.
- `filename`: sdf file (all files with extensions .vtk are accepted);
- `radius`: cylinder radius;
- `inside`: true if the sdf is negative inside.

Returns a structure containing the information of the created sdf. 
"""
function constructor_discrete_sdf(filename, radius, inside, T=Float64)
    
    # Read sdf from file
    npx, npy, npz, dx, dy, dz, dom, sdf = read_VTK_sdf(filename)
    if inside 
        field = reshape(sdf, (npx,npy,npz))
    else 
        field = reshape(-sdf, (npx,npy,npz))
    end 

    # Coordinates where the sdf values are taken 
    x0 = dom[1]
    y0 = dom[3]
    z0 = dom[5]

    xend = dom[2]
    yend = dom[4]
    zend = dom[6]

    x = range(x0; step=dx, stop=xend)
    y = range(y0; step=dy, stop=yend)
    z = range(z0; step=dz, stop=zend)

    # Build quadratic interpolation (quadratic should be enough in our case to get continuous gradient)
    itp = interpolate(field, BSpline(Quadratic(Reflect(OnCell()))))

    # Scale the interpolation on the defined coordinate grid
    sitp = scale(itp, x, y, z)  

    return SDF_Discrete{T, typeof(sitp)}(radius, sitp, dom, dx, dy, dz)  

end 

#----------------------------------
# GET INFO OF ANALYTICAL SDFs
#----------------------------------

# Get gap, gradient and hession of z-normal plane analytical SDF at point 
@inline function get_SDF_at_P_analitycal_plane_z(point::AbstractVector{T}, sdf) where T
    
    gₙ = point[3] - sdf.z0 - sdf.r
    ∂gₙ∂x = Vec3{T}(0, 0, 1)
    ∂²gₙ∂x² = zeros(Mat33{T})
    
    return gₙ, ∂gₙ∂x, ∂²gₙ∂x²
    
end 

# Get gap, gradient and hession of y-normal plane analytical SDF at point 
@inline function get_SDF_at_P_analitycal_plane_y(point::AbstractVector{T}, sdf) where T
    
    gₙ = point[2] - sdf.y0 - sdf.r
    ∂gₙ∂x = Vec3{T}(0, 1, 0)
    ∂²gₙ∂x² = zeros(Mat33{T})
    
    return -gₙ, ∂gₙ∂x, ∂²gₙ∂x²
    
end

# Get gap, gradient and hession of sphere analytical SDF at point 
@inline function get_SDF_at_P_analitycal_sphere(point::AbstractVector{T}, sdf) where T
    
    aux = Vec3(point[1] - sdf.x0, point[2] - sdf.y0, point[3] - sdf.z0)
    
    norm_aux = norm(aux)
    invnorm = 1/norm_aux
    
    gₙ = norm_aux - sdf.R - sdf.r
    ∂gₙ∂x = invnorm * aux
    ∂²gₙ∂x² = invnorm*ID3 - (invnorm^3)*(aux*aux')
    
    return gₙ, ∂gₙ∂x, ∂²gₙ∂x²
    
end

# Get gap, gradient and hession of cylinder analytical SDF at point 
@inline function get_SDF_at_P_analitycal_cylinder(point::AbstractVector{T}, sdf) where T
    
    aux = Vec3(point[1], point[2], 0)
    
    norm_aux = norm(aux)
    invnorm = 1/norm_aux
    
    gₙ = norm_aux - sdf.R + sdf.r
    ∂gₙ∂x = invnorm * aux
    ∂²gₙ∂x² = invnorm*ID3 + (invnorm^3)*(aux*aux')
    
    return -gₙ, -∂gₙ∂x, -∂²gₙ∂x²
    
end 

#----------------------------------
# INTERPOLATE DISCRETE SDF 
#----------------------------------

@inline function isinside(point, dom)

    l_x = point[1] - dom[1]  
    flag_lx = l_x>=0 && l_x<=(dom[2]-dom[1]) 
    l_y = point[2] - dom[3]
    flag_ly = l_y>=0 && l_y<=(dom[4]-dom[3])
    l_z = point[3] - dom[5] 
    flag_lz = l_z>=0 && l_z<=(dom[6]-dom[5])

    return flag_lx && flag_ly && flag_lz   

end

@inline function get_SDF_at_P_discrete(point::AbstractVector{T}, sdf) where T

    sitp = sdf.sitp

    if isinside(point, sdf.dom)

        gₙ = sitp(point...)
        ∂gₙ∂x = Interpolations.gradient(sitp, point...)
        ∂²gₙ∂x² = Interpolations.hessian(sitp, point...)

    else 

        gₙ = zero(T)
        ∂gₙ∂x = zeros(Vec3{T})
        ∂²gₙ∂x² = zeros(Mat33{T})
                
    end 

    return gₙ - sdf.r, ∂gₙ∂x, ∂²gₙ∂x²
        
end 

#----------------------------------
# CONTACT FUNCTIONS
#----------------------------------

# Get contact at point point
@inline contact_gap_specialize(point, sdf::SDF_Plane_z) = get_SDF_at_P_analitycal_plane_z(point, sdf)
@inline contact_gap_specialize(point, sdf::SDF_Sphere) = get_SDF_at_P_analitycal_sphere(point, sdf)
@inline contact_gap_specialize(point, sdf::SDF_Cylinder) = get_SDF_at_P_analitycal_cylinder(point, sdf)
@inline contact_gap_specialize(point, sdf::SDF_Plane_y) = get_SDF_at_P_analitycal_plane_y(point, sdf)
@inline contact_gap_specialize(point, sdf::SDF_Discrete) = get_SDF_at_P_discrete(point, sdf)
    
@inline function contact_gap(point, εᶜ, sdf)
        
    gₙ, ∂gₙ∂x, ∂²gₙ∂x² = contact_gap_specialize(point, sdf)
    
    pₙ, p′ₙ, Πₑ = quadratically_regularized_penalty(gₙ, εᶜ, sdf.r)
    
    return pₙ, p′ₙ, Πₑ, gₙ, ∂gₙ∂x, ∂²gₙ∂x²
    
end 



# Quadratically regulise penalty
@inline function quadratically_regularized_penalty(gₙ::T, εᶜ, r) where T
    
    ḡₙ = r/2 
    p̄ₙ = εᶜ*ḡₙ/2
    
    pₙ = zero(T)
    p′ₙ = zero(T)
    Πₑ = zero(T)
    
    # eq[111] in [2]
    if gₙ≤0

        pₙ = p̄ₙ - εᶜ*gₙ
        p′ₙ = -εᶜ
        Πₑ = (εᶜ/2)*gₙ^2 - p̄ₙ*gₙ + (εᶜ*ḡₙ^2)/6
        
    elseif gₙ≤ḡₙ
        
        aux = (εᶜ*ḡₙ-p̄ₙ)/(ḡₙ^2)
        pₙ =  aux*gₙ^2 - εᶜ*gₙ + p̄ₙ
        p′ₙ = 2*aux*gₙ - εᶜ
        Πₑ = -(εᶜ*ḡₙ-p̄ₙ)/(3*ḡₙ^2)*gₙ^3 + (εᶜ/2)*gₙ^2 - p̄ₙ*gₙ + (εᶜ*ḡₙ^2)/6
        
    end
    
    return pₙ, p′ₙ, Πₑ
    
end 


