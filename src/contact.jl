
#----------------------------------
# STRUCTURES
#----------------------------------

"""
Constructor of the structure containing the properties  of the analytical SDF of a z-normal plane.
"""
struct SDF_Plane_z{T}
    
    r::T  # beams radius
    z0::T # plane position 
    
end 

"""
Constructor of the structure containing the properties  of the analytical SDF of a y-normal plane.
"""
struct SDF_Plane_y{T}
    
    r::T  # beams radius
    y0::T # plane position 
    
end 

"""
Constructor of the structure containing the properties  of the analytical SDF of a sphere.
"""
struct SDF_Sphere{T}
    
    r::T # beams radius
    R::T # sphere radius
    x0::T # sphere position
    y0::T # sphere position
    z0::T # sphere position
    
end 

"""
Constructor of the structure containing the properties  of the analytical SDF of a z-oriented cylinder.
"""
struct SDF_Cylinder{T}
    
    r::T # beams radius
    R::T # cylinder radius
    
    
end  

# Properties of a discrete SDF
struct SDF_Discrete{T, F}
    
    r::T
    sitp::F
    dom::Array{T,1}
    dx::T
    dy::T
    dz::T
    g::T
    dg::Array{T,1}
    ddg::Array{T,2}
    
end 

#----------------------------------
# CONSTRUCTOR DISCRETE SDF
#----------------------------------

"""
Constructor of the discrete SDF from a vtk file.
"""
function constructor_discrete_sdf(filename, rWireSection, inside,  T=Float64)
    
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

    # initialise sdf variables
    g = 0
    dg = zeros(3)
    ddg = zeros(3,3)

    return SDF_Discrete{T, typeof(sitp)}(rWireSection, sitp, dom, dx, dy, dz, g, dg, ddg)  

end 

#----------------------------------
# GET INFO OF ANALYTICAL SDFs
#----------------------------------

# Get gap, gradient and hession of z-normal plane analytical SDF at P 
function get_SDF_at_P_analitycal_plane_z(P, sdf,  T=Float64)
    
    g_G = P[3] - sdf.z0 - sdf.r
    dg_G = Vec3{T}(0, 0, 1)
    ddg_G = zeros(Mat33{T})
    
    return g_G, dg_G, ddg_G
    
end 

# Get gap, gradient and hession of y-normal plane analytical SDF at P 
function get_SDF_at_P_analitycal_plane_y(P, sdf,  T=Float64)
    
    g_G = P[2] - sdf.y0 - sdf.r
    dg_G = Vec3{T}(0, 1, 0)
    ddg_G = zeros(Mat33{T})
    
    return -g_G, dg_G, ddg_G
    
end

# Get gap, gradient and hession of sphere analytical SDF at P 
function get_SDF_at_P_analitycal_sphere(P, sdf,  T=Float64)
    
    aux = Vec3{T}(P[1] .- sdf.x0, P[2] .- sdf.y0, P[3] .- sdf.z0)
    
    norm_aux = norm(aux)
    invnorm = 1/norm_aux
    
    g_G = norm_aux - sdf.R - sdf.r
    dg_G = invnorm * aux
    ddg_G = invnorm*Mat33{T}(1, 0, 0, 0, 1, 0, 0, 0, 1) + (invnorm^3)*(aux*aux')
    
    return g_G, dg_G, ddg_G
    
end

# Get gap, gradient and hession of cylinder analytical SDF at P 
function get_SDF_at_P_analitycal_cylinder(P, sdf,  T=Float64)
    
    aux = Vec3{T}(P[1], P[2], 0)
    
    norm_aux = norm(aux)
    invnorm = 1/norm_aux
    
    g_G = norm_aux - sdf.R + sdf.r
    dg_G = invnorm * aux
    ddg_G = invnorm*I + (invnorm^3)*(aux*aux')
    
    return -g_G, -dg_G, -ddg_G
    
end 

#----------------------------------
# INTERPOLATE DISCRETE SDF 
#----------------------------------

function get_SDF_at_P_discrete(xP, sdf,  T=Float64)
    
    dom = sdf.dom
    sitp = sdf.sitp

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

    return g - sdf.r, dg, ddg
        
end 

#----------------------------------
# CONTACT FUNCTIONS
#----------------------------------

# Get contact at point GP
function get_contact_GP(GP, epsC, sdf, T=Float64)
    
    if typeof(sdf) == SDF_Plane_z{T}
        
        g_G, dg_G, ddg_G = get_SDF_at_P_analitycal_plane_z(GP, sdf, T)
        
    elseif typeof(sdf) == SDF_Sphere{T}
        
        g_G, dg_G, ddg_G = get_SDF_at_P_analitycal_sphere(GP, sdf, T)
        
    elseif typeof(sdf) == SDF_Cylinder{T}
        
        g_G, dg_G, ddg_G = get_SDF_at_P_analitycal_cylinder(GP, sdf, T)
        
    elseif typeof(sdf) == SDF_Plane_y{T}
        
        g_G, dg_G, ddg_G = get_SDF_at_P_analitycal_plane_y(GP, sdf, T)
        
    else
        g_G, dg_G, ddg_G = get_SDF_at_P_discrete(GP, sdf, T)
    end 
    
    fc_eps, dfc_eps, Pic_eps = quadratically_regularized_penalty(g_G, epsC, sdf.r, T)
    
    return fc_eps, dfc_eps, Pic_eps, g_G, dg_G, ddg_G
    
end 

# Quadratically regulise penalty
function quadratically_regularized_penalty(g, epsC, r,  T=Float64)
    
    gbar = 0.1*r
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


