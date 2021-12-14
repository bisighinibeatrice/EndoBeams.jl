
#----------------------------------
# STRUCTURES
#----------------------------------

"""
    SDF = struct SDF_Plane_z{T}

Constructor of the structure containing the properties of the analytical SDF of a xy plane:
- `r`: beam radius (`::T`);
- `z0`: plane position along z (`::T`).

Returns a structure containing the information of the created SDF. 
"""
struct SDF_Plane_z{T}
    
    r::T  # beams radius
    z0::T # plane position 
    
end 

"""
    SDF = struct SDF_Plane_y{T}

Constructor of the structure containing the properties  of the analytical SDF of a xz plane:
- `r`: beam radius (`::T`);
- `y0`: plane position along y (`::T`).

Returns a structure containing the information of the created SDF.  
"""
struct SDF_Plane_y{T}
    
    r::T  # beams radius
    y0::T # plane position 
    
end 

"""
    SDF = struct SDF_Sphere{T}

Constructor of the structure containing the properties  of the analytical SDF of a sphere:
- `r`: beam radius (`::T`);
- `R`: sphere radius (`::T`);
- `x0`: sphere centre position along x (`::T`);
- `y0`: sphere centre position along y (`::T`);
- `z0`: sphere centre position along z (`::T`).
    
Returns a structure containing the information of the created SDF. 
"""
struct SDF_Sphere{T}
    
    r::T # beams radius
    R::T # sphere radius
    x0::T # sphere position
    y0::T # sphere position
    z0::T # sphere position
    
end 

"""
    SDF = struct SDF_Cylinder{T}

Constructor of the structure containing the properties of the analytical SDF of an infinite cylinder oriented along z:
- `r`: beam radius (`::T`);
- `R`: cylinder radius (`::T`).
    
Returns a structure containing the information of the created SDF. 
"""
struct SDF_Cylinder{T}
    
    r::T # beams radius
    R::T # cylinder radius
    
    
end  

# Properties of a discrete SDF
struct SDF_Discrete{T, F}
    
    r::T
    sitp::F
    dom::Vector{T}
    dx::T
    dy::T
    dz::T
    g::T
    dg::Vector{T}
    ddg::Matrix{T}
    
end 

#----------------------------------
# CONSTRUCTOR DISCRETE SDF
#----------------------------------
"""
    SDF = constructor_discrete_SDF(filename, rWireSection, inside,  T=Float64)

Constructor of the discrete SDF from a vtk file:
- `filename`: SDF file (files with extensions .vtk are accepted);
- `rWireSection`: cylinder radius (`::T`);
- `inside`: true if the SDF is negative inside (`::Bool`).

Returns a structure containing the information relative to the created SDF. 
"""
function constructor_discrete_SDF(filename, rWireSection, inside,  T=Float64)
    
    # Read SDF from file
    npx, npy, npz, dx, dy, dz, dom, SDF = read_VTK_sdf(filename)
    if inside 
        field = reshape(SDF, (npx,npy,npz))
    else 
        field = reshape(-SDF, (npx,npy,npz))
    end 

    # Coordinates where the SDF values are taken 
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

    # initialise SDF variables
    g = 0
    dg = zeros(3)
    ddg = zeros(3,3)

    return SDF_Discrete{T, typeof(sitp)}(rWireSection, sitp, dom, dx, dy, dz, g, dg, ddg)  

end 

#----------------------------------
# GET INFO OF ANALYTICAL SDFs
#----------------------------------

# Get gap, gradient and hession of z-normal plane analytical SDF at P 
function get_SDF_at_P_analitycal_plane_z(P, SDF,  T=Float64)
    
    g_G = P[3] - SDF.z0 - SDF.r
    dg_G = Vec3(0, 0, 1)
    ddg_G = zeros(Mat33{T})
    
    return g_G, dg_G, ddg_G
    
end 

# Get gap, gradient and hession of y-normal plane analytical SDF at P 
function get_SDF_at_P_analitycal_plane_y(P, SDF,  T=Float64)
    
    g_G = P[2] - SDF.y0 - SDF.r
    dg_G = Vec3(0, 1, 0)
    ddg_G = zeros(Mat33{T})
    
    return -g_G, dg_G, ddg_G
    
end

# Get gap, gradient and hession of sphere analytical SDF at P 
function get_SDF_at_P_analitycal_sphere(P, SDF,  T=Float64)
    
    aux = Vec3(P[1] - SDF.x0, P[2] - SDF.y0, P[3] - SDF.z0)
    
    norm_aux = norm(aux)
    invnorm = 1/norm_aux
    
    g_G = norm_aux - SDF.R - SDF.r
    dg_G = invnorm * aux
    ddg_G = invnorm*ID3 + (invnorm^3)*(aux*aux')
    
    return g_G, dg_G, ddg_G
    
end

# Get gap, gradient and hession of cylinder analytical SDF at P 
function get_SDF_at_P_analitycal_cylinder(P, SDF,  T=Float64)
    
    aux = Vec3(P[1], P[2], 0)
    
    norm_aux = norm(aux)
    invnorm = 1/norm_aux
    
    g_G = norm_aux - SDF.R + SDF.r
    dg_G = invnorm * aux
    ddg_G = invnorm*I + (invnorm^3)*(aux*aux')
    
    return -g_G, -dg_G, -ddg_G
    
end 

#----------------------------------
# INTERPOLATE DISCRETE SDF 
#----------------------------------

function get_SDF_at_P_discrete(xP, SDF,  T=Float64)
    
    dom = SDF.dom
    sitp = SDF.sitp

    l_x = xP[1] - dom[1]  
    flag_lx = l_x>=0 && l_x<=(dom[2]-dom[1]) 
    l_y = xP[2] - dom[3]
    flag_ly = l_y>=0 && l_y<=(dom[4]-dom[3])
    l_z = xP[3] - dom[5] 
    flag_lz = l_z>=0 && l_z<=(dom[6]-dom[5])

    flag_in = flag_lx*flag_ly*flag_lz               

    
    if flag_in == 0
        
        g = zero(T)
        dg = zeros(Vec3{T})
        ddg = zeros(Mat33{T})

    else 
    
        g = sitp(xP...)
        dg = Interpolations.gradient(sitp, xP...)
        ddg = Interpolations.hessian(sitp, xP...)
                
    end 

    return g - SDF.r, dg, ddg
        
end 

#----------------------------------
# CONTACT FUNCTIONS
#----------------------------------

# Get contact at point GP
get_contact_GP_specialize(GP, SDF::SDF_Plane_z, T) = get_SDF_at_P_analitycal_plane_z(GP, SDF, T)
get_contact_GP_specialize(GP, SDF::SDF_Sphere, T) = get_SDF_at_P_analitycal_sphere(GP, SDF, T)
get_contact_GP_specialize(GP, SDF::SDF_Cylinder, T) = get_SDF_at_P_analitycal_cylinder(GP, SDF, T)
get_contact_GP_specialize(GP, SDF::SDF_Plane_y, T) = get_SDF_at_P_analitycal_plane_y(GP, SDF, T)
get_contact_GP_specialize(GP, SDF::SDF_Discrete, T) = get_SDF_at_P_discrete(GP, SDF, T)
    
function get_contact_GP(GP, epsC, SDF, T=Float64)
        
    g_G, dg_G, ddg_G = get_contact_GP_specialize(GP, SDF, T)
    
    fc_eps, dfc_eps, Pic_eps = quadratically_regularized_penalty(g_G, epsC, SDF.r, T)
    
    return fc_eps, dfc_eps, Pic_eps, g_G, dg_G, ddg_G
    
end 

# Quadratically regulise penalty
function quadratically_regularized_penalty(g, epsC, r,  T=Float64)
    
    gbar = r/10
    fbar = epsC*gbar/2
    aux = (epsC*gbar-fbar)/(gbar^2)
     
    fc_eps = zero(T)
    dfc_eps = zero(T)
    Pic_eps = zero(T)
    
    # eq[111] in [2]
    if g<=0

        fc_eps = fbar - epsC*g
        dfc_eps = -epsC
        Pic_eps = (epsC/2)*g^2 - fbar*g + (epsC*gbar^2)/6
        
    elseif gbar>=g && g>0
        
        fc_eps =  aux*g^2 - epsC*g + fbar
        dfc_eps = 2*aux*g - epsC
        Pic_eps = -(epsC*gbar-fbar)/(3*gbar^2)*g^3 + (epsC/2)*g^2 - fbar*g + (epsC*gbar^2)/6
        
    end
    
    return fc_eps, dfc_eps, Pic_eps
    
end 


