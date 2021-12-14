#----------------------------------
# STRUCTURE
#----------------------------------

struct MyBeam{T}
    
    ind::Int # index of this beam 
    indGP::Vec3{Int} #indices of the GP associated to this beam
    node1::Int # index of node 1   
    node2::Int # index of node 2
    l0::T # initial beam length
    R0::Mat33{T} # initial beam rotation matrix
    Kint::Mat77{T} # beam internal matrix
    numberInterpolationPoints::Int # number of interpolation points (-> visualization)
    sparsity_map:: Array{Int,1} # sparsity map from local indices (beam matrix) to global indices (gloabl sparse matrix) -> computed in the constructor of the sparse matrices

end

#----------------------------------
# CONSTRUCTORS
#----------------------------------
"""
    beams = constructor_beams(allnodes, connectivity, mat, geom, numberInterpolationPoints, Re0=nothing, T=Float64)

Constructor of the beams StructArray:
- `allnodes`: nodes (`::StructArray{MyNode}` created with 'constructor_nodes');
- `connectivity`: mesh connectivity (`::Vec2{Int}`);
- `mat`: mesh material properties (`::Material{T}`);
- `geom`: mesh geometrical properties (`::Geometry{T}`);
- `numberInterpolationPoints`: number of points used for the interpolation of the beam centerline (`::Int`);
- `Re0`: (not mandatory) initial rotation matrix of the beam elements (`::Mat33`).

Returns a `StructArray{MyBeam}` containing the information relative to the beam elements. 
"""
function constructor_beams(allnodes, connectivity, mat, geom, numberInterpolationPoints, Re0=nothing, T=Float64)
    
    if isnothing(Re0)
        beams = StructArray(constructor_beam(i, (i-1)*3 .+ [1,2,3], allnodes[connectivity[i][1]], allnodes[connectivity[i][2]], mat, geom, numberInterpolationPoints, T) for i in 1:size(connectivity,1))
    else 
        beams = StructArray(constructor_beam_Re0(i, (i-1)*3 .+ [1,2,3], allnodes[connectivity[i][1]], allnodes[connectivity[i][2]], mat, geom, numberInterpolationPoints, Re0[i], T) for i in 1:size(connectivity,1))
    end 

    return beams
    
end 

# Constructor of the one beam (MyBeam)
function constructor_beam(ind, indGP, node1, node2, mat, geom, numberInterpolationPoints, T=Float64)
    
    i1 = node1.i
    i2 = node2.i
    R0 = get_R0_beam(node1, node2, T)
    l0 = norm(node1.pos - node2.pos)
    Kint = get_Kint_beam(mat, geom, l0, T)
    
    return MyBeam{T}(ind, indGP, i1, i2, l0, R0, Kint, numberInterpolationPoints, zeros(Int, 144))
    
end 

# Constructor of the one beam (MyBeam) given initial rotation
function constructor_beam_Re0(ind, indGP, node1, node2, mat, geom, numberInterpolationPoints, Re0, T=Float64)
    
    i1 = node1.i   
    i2 = node2.i  
    R0 = Re0
    l0 = norm(node1.pos - node2.pos)   
    Kint = get_Kint_beam(mat, geom, l0, T)
    
    return MyBeam{T}(ind, indGP, i1, i2, l0, R0, Kint, numberInterpolationPoints, zeros(Int, 144))
    
end 

#----------------------------------
# UTILS 
#----------------------------------

# Compute R0 matrix
function get_R0_beam(node1, node2, T=Float64)
    
    tol =  2*eps(T)
    
    X1 = node1.pos
    X2 = node2.pos
    
    l0 = norm(X2-X1)
    
    E1 = Vec3(1, 0, 0)
    
    E1_0 = (X2-X1)/l0
    
    v = cross(E1, E1_0)
    
    s = norm(v)
    c = dot(E1, E1_0)
    
    if (c<-(1-tol))
        
        R0 = Mat33{T}(-1, 0, 0, 0, 1, 0, 0, 0, -1)
        
    elseif (c > (1-tol))
        
        R0 = Mat33{T}(1, 0, 0, 0, 1, 0, 0, 0, 1)
        
    else 
        
        Sv = get_skew_skymmetric_matrix_from_vector(v)
        R0 = ID3 + Sv + ((1-c)/s^2)*Sv*Sv
    end 
    
    return R0
    
end 

#  Compute Kint matrix
function get_Kint_beam(mat, geom, l0, T=Float64)
    
    Kint_bar = Mat77{T}(geom.A*mat.E/l0, 0, 0, 0, 0, 0, 0,
    0, mat.G*geom.J/l0, 0, 0, -mat.G*geom.J/l0, 0, 0,
    0, 0, 4*mat.E*geom.I33/l0, 0, 0, 2*mat.E*geom.I33/l0, 0,
    0, 0, 0, 4*mat.E*geom.I22/l0, 0, 0, 2*mat.E*geom.I22/l0,
    0, -mat.G*geom.J/l0, 0, 0, mat.G*geom.J/l0, 0, 0,
    0, 0, 2*mat.E*geom.I33/l0, 0, 0, 4*mat.E*geom.I33/l0, 0,  
    0, 0, 0, 2*mat.E*geom.I22/l0, 0, 0, 4*mat.E*geom.I22/l0)
    
    return Kint_bar
    
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
        Re0 = e.R0
        
        # information from node 1 and 2
        X1, X2 = get_local_pos(e, allnodes)
        uk1, uk2 = get_local_displ(e, allnodes)
        R1, R2 = get_local_rot(e, allnodes)
        
        # current position of node i1 (global rs)
        x1 =  X1 + uk1
        x2 =  X2 + uk2
        
        # -------------------------------------------------------------------------------------------
        # rigidly_moving rs of the deformed configuration (v1, v2, v3)
        # -------------------------------------------------------------------------------------------
        
        l0 = e.l0        
        ln = norm(x2-x1)
        
        v1 = (x2 - x1) / ln
        
        t2_1 = R1 * Re0 * Vec3(0, 1, 0)
        t2_2 = R2 * Re0 * Vec3(0, 1, 0)  
        p = (t2_1 + t2_2)/2
        
        v3 = cross(v1, p)
        v3 = v3 / norm(v3)
        v2 = cross(v3, v1)
        
        Re = [v1 v2 v3]
        
        R1_bar = Re' * R1 * Re0
        R2_bar = Re' * R2 * Re0
        
        psil1 = get_angle_from_rotation_matrix(R1_bar)
        psil2 = get_angle_from_rotation_matrix(R2_bar)

        N1(xi) = 1-xi/l0
        N2(xi) = 1-N1(xi)
        N3(xi) = xi.*(1-xi/l0).^2
        N4(xi) = -(1-xi/l0).*(xi.^2/l0)
        
        dxi = l0/allbeams.numberInterpolationPoints[i]
        xi_list = 0:dxi:l0
        if length(xi_list) != allbeams.numberInterpolationPoints[i]+1
            xi_list = range(0,stop=l0,length=allbeams.numberInterpolationPoints[i]+1)
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

            xOG_j = N1_i*x1 + N2_i*x2 + Re*ul

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
