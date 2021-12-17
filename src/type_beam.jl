#----------------------------------
# STRUCTURE
#----------------------------------

struct MyBeam{T}
    
    ind::Int # index of this beam 
    indGP::Vec3{Int} #indices of the GP associated to this beam
    node1::Int # index of node 1   
    node2::Int # index of node 2
    l₀::T # initial beam length
    R₀::Mat33{T} # initial beam rotation matrix
    Kint::Mat77{T} # beam internal matrix
    numberInterpolationPoints::Int # number of interpolation points (-> visualization)
    sparsity_map:: Array{Int,1} # sparsity map from local indices (beam matrix) to global indices (gloabl sparse matrix) -> computed in the constructor of the sparse matrices

end

#----------------------------------
# CONSTRUCTORS
#----------------------------------
"""
    beams = constructor_beams(allnodes, connectivity, mat, geom, numberInterpolationPoints, Rₑ⁰=nothing, T=Float64)

Constructor of the beams StructArray:
- `allnodes`: nodes StructArray (created with constructor_nodes);
- `connectivity`: connectivity of the mesh (Vec2{Int});
- `mat`: struct containing the material properties of the mesh (Material{T});
- `geom`: struct containing the geomtrical properties of the mesh (Geometry{T});
- `numberInterpolationPoints`: pnumber of points used for the interpolation of the beam centreline (Int);
- `Rₑ⁰`: (not mandatory) initial rotation of the beam elements.

Returns a StructArray{MyBeam}, structure containing the information of the beam elements. 
"""
function constructor_beams(allnodes, connectivity, mat, geom, numberInterpolationPoints, Rₑ⁰=nothing, T=Float64)
    
    if isnothing(Rₑ⁰)
        beams = StructArray(constructor_beam(i, (i-1)*3 .+ [1,2,3], allnodes[connectivity[i][1]], allnodes[connectivity[i][2]], mat, geom, numberInterpolationPoints, T) for i in 1:size(connectivity,1))
    else 
        beams = StructArray(constructor_beam_Re0(i, (i-1)*3 .+ [1,2,3], allnodes[connectivity[i][1]], allnodes[connectivity[i][2]], mat, geom, numberInterpolationPoints, Rₑ⁰[i], T) for i in 1:size(connectivity,1))
    end 

    return beams
    
end 

# Constructor of the one beam (MyBeam)
function constructor_beam(ind, indGP, node1, node2, mat, geom, numberInterpolationPoints, T=Float64)
    
    i1 = node1.i
    i2 = node2.i
    R₀ = get_R0_beam(node1, node2, T)
    l₀ = norm(node1.pos - node2.pos)
    Kint = K̄ᵢₙₜ_beam(mat, geom, l₀, T)
    
    return MyBeam{T}(ind, indGP, i1, i2, l₀, R₀, Kint, numberInterpolationPoints, zeros(Int, 144))
    
end 

# Constructor of the one beam (MyBeam) given initial rotation
function constructor_beam_Re0(ind, indGP, node1, node2, mat, geom, numberInterpolationPoints, Rₑ⁰, T=Float64)
    
    i1 = node1.i   
    i2 = node2.i  
    R₀ = Rₑ⁰
    l₀ = norm(node1.pos - node2.pos)   
    Kint = K̄ᵢₙₜ_beam(mat, geom, l₀, T)
    
    return MyBeam{T}(ind, indGP, i1, i2, l₀, R₀, Kint, numberInterpolationPoints, zeros(Int, 144))
    
end 

#----------------------------------
# UTILS 
#----------------------------------

# Compute R₀ matrix
function get_R0_beam(node1, node2, T=Float64)
    
    tol =  2*eps(T)
    
    X₁ = node1.pos
    X₂ = node2.pos
    
    l₀ = norm(X₂-X₁)
    
    E1 = Vec3(1, 0, 0)
    
    E1_0 = (X₂-X₁)/l₀
    
    v = cross(E1, E1_0)
    
    s = norm(v)
    c = dot(E1, E1_0)
    
    if (c<-(1-tol))
        
        R₀ = Mat33{T}(-1, 0, 0, 0, 1, 0, 0, 0, -1)
        
    elseif (c > (1-tol))
        
        R₀ = Mat33{T}(1, 0, 0, 0, 1, 0, 0, 0, 1)
        
    else 
        
        Sv = skew_skymmetric_matrix_from_vector(v)
        R₀ = ID3 + Sv + ((1-c)/s^2)*Sv*Sv
    end 
    
    return R₀
    
end 


#  Compute Kint matrix
function K̄ᵢₙₜ_beam(mat, geom, l₀)
    
    K̄ᵢₙₜū = geom.A*mat.E/l₀
    K̄ᵢₙₜΘ̅ᵢᵢ = @SVector [mat.G*geom.J/l₀, 4*mat.E*geom.I33/l₀, 4*mat.E*geom.I22/l₀]
    K̄ᵢₙₜΘ̅ᵢⱼ = @SVector [-mat.G*geom.J/l₀, 2*mat.E*geom.I33/l₀, 2*mat.E*geom.I22/l₀]
    
    return K̄ᵢₙₜū, K̄ᵢₙₜΘ̅ᵢᵢ, K̄ᵢₙₜΘ̅ᵢⱼ
    
end

# Obtain beam centerline position by interpolation
function get_centerline!(positions, connectivity, allnodes, allbeams)
    
    Xcount = 0
    index = 0

    for (i, e) in enumerate(allbeams)
        
        # -------------------------------------------------------------------------------------------
        # information from node 1 and 2
        # -------------------------------------------------------------------------------------------
        # retrieve the matrix Re_0 of the beam
        Rₑ⁰ = e.R₀
        
        # information from node 1 and 2
        X₁, X₂ = local_pos(e, allnodes)
        u₁, u₂ = local_disp(e, allnodes)
        R₁, R₂ = local_rot(e, allnodes)
        
        # current position of node i1 (global rs)
        x₁ =  X₁ + u₁
        x₂ =  X₂ + u₂
        
        # -------------------------------------------------------------------------------------------
        # rigidly_moving rs of the deformed configuration (v1, v2, v3)
        # -------------------------------------------------------------------------------------------
        
        l₀ = e.l₀        
        lₙ = norm(x₂-x₁)
        
        v1 = (x₂ - x₁) / lₙ
        
        t2_1 = R₁ * Rₑ⁰ * Vec3(0, 1, 0)
        t2_2 = R₂ * Rₑ⁰ * Vec3(0, 1, 0)  
        p = (t2_1 + t2_2)/2
        
        v3 = cross(v1, p)
        v3 = v3 / norm(v3)
        v2 = cross(v3, v1)
        
        Rₑ = [v1 v2 v3]
        
        R̅₁ = Rₑ' * R₁ * Rₑ⁰
        R̅₂ = Rₑ' * R₂ * Rₑ⁰
        
        psil1 = angle_from_rotation_matrix(R̅₁)
        psil2 = angle_from_rotation_matrix(R̅₂)

        N1(xi) = 1-xi/l₀
        N2(xi) = 1-N1(xi)
        N3(xi) = xi.*(1-xi/l₀).^2
        N4(xi) = -(1-xi/l₀).*(xi.^2/l₀)
        
        dxi = l₀/allbeams.numberInterpolationPoints[i]
        xi_list = 0:dxi:l₀
        if length(xi_list) != allbeams.numberInterpolationPoints[i]+1
            xi_list = range(0,stop=l₀,length=allbeams.numberInterpolationPoints[i]+1)
        end 

        for x_i in xi_list
            
            index = index + 1

            N1_i = N1(x_i)
            N2_i = N2(x_i)
            N3_i = N3(x_i)
            N4_i = N4(x_i)

            P1 = Mat36(
                0, 0, 0, 
                0, 0, -N3_i, 
                0, N3_i, 0, 
                0, 0, 0,
                0, 0, -N4_i,
                0, N4_i, 0)
            
            ul = P1* Vec6(psil1[1], psil1[2], psil1[3],  psil2[1], psil2[2], psil2[3]) # local cross section displacement, eq[38] Le 2014

            xOG_j = N1_i*x₁ + N2_i*x₂ + Rₑ*ul

            positions[index, 1] = xOG_j[1]
            positions[index, 2] = xOG_j[2]
            positions[index, 3] = xOG_j[3]
            
        end 
        
        Xcount_ini = Xcount + 1
        aux =  1:length(xi_list)
        Xcount_end = Xcount + length(aux)
        Xcount = Xcount_end
        ind_nodes_elem = Xcount_ini:1:Xcount_end

        connectivity[:, i] =  ind_nodes_elem
        
    end

end
