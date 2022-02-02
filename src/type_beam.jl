#----------------------------------
# STRUCTURE
#----------------------------------

struct Beam{T, Tᵏ}
    
    ind::Int # index of this beam 
    indGP::Vec3{Int} #indices of the GP associated to this beam
    node1::Int # index of node 1   
    node2::Int # index of node 2
    l₀::T # initial beam length
    Rₑ⁰::Mat33{T} # initial beam rotation matrix
    K̄ᵢₙₜ::Tᵏ # beam internal matrix
    numberInterpolationPoints::Int # number of interpolation points (-> visualization)
    sparsity_map::Vector{Int} # sparsity map from local indices (beam matrix) to global indices (gloabl sparse matrix) -> computed in the constructor of the sparse matrices

end

#----------------------------------
# CONSTRUCTORS
#----------------------------------
"""
    beams = constructor_beams(nodes, connectivity, mat, geom, numberInterpolationPoints, Rₑ⁰=nothing, T=Float64)

Constructor of the beams StructArray:
- `nodes`: nodes StructArray (created with constructor_nodes);
- `connectivity`: connectivity of the mesh (Vec2{Int});
- `mat`: struct containing the material properties of the mesh (Material{T});
- `geom`: struct containing the geomtrical properties of the mesh (Geometry{T});
- `numberInterpolationPoints`: pnumber of points used for the interpolation of the beam centreline (Int);
- `Rₑ⁰`: (not mandatory) initial rotation of the beam elements.

Returns a StructArray{Beam}, structure containing the information of the beam elements. 
"""

function constructor_beams(nodes, connectivity, mat, geom, numberInterpolationPoints, Rₑ⁰=nothing)
    
    beams = StructArray(constructor_beam(i, (i-1)*3 .+ [1,2,3], nodes[connectivity[i][1]], nodes[connectivity[i][2]], mat, geom, numberInterpolationPoints, Rₑ⁰ isa AbstractVector ? Rₑ⁰[i] : Rₑ⁰) for i in 1:size(connectivity,1))

    return beams
    
end 

    
# Constructor of the one beam (Beam) given initial rotation
function constructor_beam(ind, indGP, node1::Node{T}, node2::Node{T}, mat, geom, numberInterpolationPoints, Rₑ⁰) where T
    
    i1 = node1.i   
    i2 = node2.i  
    l₀ = norm(node1.X₀ - node2.X₀)   
    K̄ᵢₙₜ = K̄ᵢₙₜ_beam(mat, geom, l₀)
    
    return Beam{T, typeof(K̄ᵢₙₜ)}(ind, indGP, i1, i2, l₀, Rₑ⁰, K̄ᵢₙₜ, numberInterpolationPoints, zeros(Int, 144))
    
end 

constructor_beam(ind, indGP, node1::Node{T}, node2::Node{T}, mat, geom, numberInterpolationPoints, Rₑ⁰::Nothing) where T = constructor_beam(ind, indGP, node1, node2, mat, geom, numberInterpolationPoints, local_R⁰(node1.X₀, node2.X₀))


# Obtain beam centerline position by interpolation
function get_centerline!(positions, connectivity, nodes, allbeams)
    
    Xcount = 0
    index = 0

    for e in LazyRows(allbeams)

        i = e.ind
        
        # -------------------------------------------------------------------------------------------
        # information from node 1 and 2
        # -------------------------------------------------------------------------------------------
        # retrieve the matrix Rₑ⁰ of the beam
        Rₑ⁰ = e.Rₑ⁰
        
        # information from node 1 and 2
        X₀₁, X₀₂ = nodes.X₀[e.node1], nodes.X₀[e.node2]
        u₁, u₂ = nodes.u[e.node1], nodes.u[e.node2]
        R₁, R₂ = nodes.R[e.node1], nodes.R[e.node2]
        
        # current position of node i1 (global rs)
        x₁ =  X₀₁ + u₁
        x₂ =  X₀₂ + u₂
        
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
        
        psil1 = toangle(R̅₁)
        psil2 = toangle(R̅₂)

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
