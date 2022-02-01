#-------------------------------------------------------------------
# FUNCTIONS AT THE BEGINNING OF EACH BEAM CONTRIBUTIONS COMPUTATION
#-------------------------------------------------------------------

# Gets the position of the two nodes forming the given beam element
@inline function local_pos(e, allnodes)
    
    return allnodes.pos[e.node1], allnodes.pos[e.node2]
    
end 

# Gets the displacement of the two nodes forming the given beam element
@inline function local_disp(e, allnodes)
    
    return allnodes.u[e.node1], allnodes.u[e.node2]
    
end 

# Gets the velocity of the two nodes forming the given beam element
@inline function local_vel(e, allnodes)
    
    return allnodes.udt[e.node1], allnodes.udt[e.node2]
    
end 

# Gets the acceleration of the two nodes forming the given beam element
@inline function local_acc(e, allnodes)
    
    return allnodes.udtdt[e.node1], allnodes.udtdt[e.node2]
    
end 

# Gets the angular velocity of the two nodes forming the given beam element
@inline function local_ang_vel(e, allnodes)
    
    return allnodes.wdt[e.node1], allnodes.wdt[e.node2]
    
end 

# Gets the angular acceleration of the two nodes forming the given beam element
@inline function local_ang_acc(e, allnodes)
    
    return allnodes.wdtdt[e.node1], allnodes.wdtdt[e.node2]
    
end 

# Gets the rotation matrix of the two nodes forming the given beam element
@inline function local_rot(e, allnodes)
    
    return allnodes.R[e.node1], allnodes.R[e.node2]
    
end 

# Gets the rotation matrix displacement of the two nodes forming the given beam element
@inline function local_rot_delt(e, allnodes)
    
    return allnodes.Delt[e.node1], allnodes.Delt[e.node2]
    
end 

#----------------------------------
# RESHAPE FUNCTIONS
#----------------------------------

function divide_Vec12_into_Vec4Vec3(a::Vec12)
    
    f = Vec3(a[1], a[2], a[3])
    e = Vec3(a[4], a[5], a[6])
    c = Vec3(a[7], a[8], a[9])
    d = Vec3(a[10], a[11], a[12])
    
    return Vec4(f,  e,  c,  d)
    
end

function compute_ssmatrices_Vec4Vec3_return_Vec4Mat33(a)
    
    aux1 = skew_skymmetric_matrix_from_vector(a[1])
    aux2 = skew_skymmetric_matrix_from_vector(a[2])
    aux3 = skew_skymmetric_matrix_from_vector(a[3])
    aux4 = skew_skymmetric_matrix_from_vector(a[4])
    
    return Vec4(aux1, aux2, aux3, aux4)
end

function shape_Vec12(a, b, c, d)
    
    return Vec12(a[1], a[2], a[3], b[1], b[2], b[3], c[1], c[2], c[3], d[1], d[2], d[3])
    
end

function divide_Vec12_into_Vec4Vec3(a)
    
    f = Vec3(a[1], a[2], a[3])
    e = Vec3(a[4], a[5], a[6])
    c = Vec3(a[7], a[8], a[9])
    d = Vec3(a[10], a[11], a[12])
    
    return Vec4(f,  e,  c,  d)
    
end

function reshapeA1!(M, v)
    
    M[2,1] = M[2,1] - (v[8]-v[2])
    M[3,1] = M[3,1] - (v[9]-v[3])
    
    M[2,7] = M[2,7] + (v[8]-v[2])
    M[3,7] = M[3,7] + (v[9]-v[3])
    
end

function reshapeA2!(M, v)
    
    M[2,1] = M[2,1] - (v[3]-v[9])
    M[3,1] = M[3,1] - (v[8]-v[2])
    
    M[2,7] = M[2,7] + (v[3]-v[9])
    M[3,7] = M[3,7] + (v[8]-v[2])
    
end

function reshapeP1G(N3, N4, Θ̅₁, Θ̅₂)
    
    return Vec3(0, N3*Θ̅₁[3]+N4*Θ̅₂[3], -(N3*Θ̅₁[2]+N4*Θ̅₂[2]))
    
end 

#------------------------------------------
# FUNCTIONS TO COMPUTE SPECIFIC MATRICES
#------------------------------------------



function compute_E_or_Et(Rₑ)
    
    Re11 = Rₑ[1,1]
    Re21 = Rₑ[2,1]
    Re31 = Rₑ[3,1]
    
    Re12 = Rₑ[1,2]
    Re22 = Rₑ[2,2]
    Re32 = Rₑ[3,2]
    
    Re13 = Rₑ[1,3]
    Re23 = Rₑ[2,3]
    Re33 = Rₑ[3,3]
    
    E = Mat1212(
    Re11, Re21, Re31, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    Re12, Re22, Re32, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    Re13, Re23, Re33, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, Re11, Re21, Re31, 0, 0, 0, 0, 0, 0,
    0, 0, 0, Re12, Re22, Re32, 0, 0, 0, 0, 0, 0,
    0, 0, 0, Re13, Re23, Re33, 0, 0, 0, 0, 0, 0,
    0, 0, 0,  0, 0, 0, Re11, Re21, Re31, 0, 0, 0,
    0, 0, 0,  0, 0, 0, Re12, Re22, Re32, 0, 0, 0,
    0, 0, 0,  0, 0, 0, Re13, Re23, Re33, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, Re11, Re21, Re31,
    0, 0, 0, 0, 0, 0, 0, 0, 0, Re12, Re22, Re32,
    0, 0, 0, 0, 0, 0, 0, 0, 0, Re13, Re23, Re33)
    
    return E
end

function compute_r(v)
    
    r = Vec12(-v[1], -v[2], -v[3], 0, 0, 0, v[1], v[2], v[3], 0, 0, 0)
    
    return r
    
end

function compute_GT_P_eta(lₙ, R̅₁, R̅₂)
    
    # eqC7 in [2]: η, eta11, eta12, eta21, eta22
    t1 = R̅₁[:,2]
    tl = R̅₂[:,2]
    p = (v1+v2)/2
    q1 = p[1]; q2 = p[2]
    q11 = t1[1]; q12 = t1[2]
    q21 = tl[1]; q22 = tl[2]
    
    η = q1/q2
    eta11 = q11/q2
    eta12 = q12/q2
    eta21 = q21/q2
    eta22 = q22/q2
    
    # eqC6 in [2]: GT
    GTu1 = @SMatrix [0 0     η/lₙ;
                     0 0     1/lₙ  ;
                     0 -1/lₙ 0      ]
    
    GTΘ1 = @SMatrix [eta12/2 -eta11/2 0;
                     0       0        0;
                     0       0        0]

    GTu2 = -GTu1

    GTΘ1 = @SMatrix [eta22/2 -eta21/2 0;
                     0       0        0;
                     0       0        0]

    Pu1 = [-GTu1; -GTu1]

    PΘ1 = [ID3-GTΘ1; -GTΘ1]

    Pu2 = [-GTu2; -GTu2]

    PΘ2 = [-GTΘ1; ID3-GTΘ1]
    
    return GTu1, GTΘ1, GTu2, GTΘ2, Pu1, PΘ1, Pu2, PΘ2, η
    
end

function compute_Kh_bar(Kh1_bar, Kh2_bar) 
    
    Kh_bar = Mat77(
    0, 0, 0, 0, 0, 0, 0,
    0, Kh1_bar[1,1], Kh1_bar[2,1], Kh1_bar[3,1], 0, 0, 0,
    0, Kh1_bar[1,2], Kh1_bar[2,2], Kh1_bar[3,2], 0, 0, 0,
    0, Kh1_bar[1,3], Kh1_bar[2,3], Kh1_bar[3,3], 0, 0, 0,
    0, 0, 0, 0, Kh2_bar[1,1], Kh2_bar[2,1], Kh1_bar[3,1],
    0, 0, 0, 0, Kh2_bar[1,2], Kh2_bar[2,2], Kh2_bar[3,2],
    0, 0, 0, 0, Kh2_bar[1,3], Kh2_bar[2,3], Kh2_bar[3,3])  
    
    return Kh_bar
    
end

function compute_B_bar(Tₛ⁻¹Θ̅₁, Tₛ⁻¹Θ̅₂) 
    
    return Tₛ⁻¹Θ̅₁, Tₛ⁻¹Θ̅₂
    
end

function compute_B_star_bar(rT, PET) 
    
    B_star_bar = Mat712(
    rT[1], PET[1,1], PET[2,1], PET[3,1], PET[4,1], PET[5,1], PET[6,1],
    rT[2], PET[1,2], PET[2,2], PET[3,2], PET[4,2], PET[5,2], PET[6,2],
    rT[3], PET[1,3], PET[2,3], PET[3,3], PET[4,3], PET[5,3], PET[6,3],
    rT[4], PET[1,4], PET[2,4], PET[3,4], PET[4,4], PET[5,4], PET[6,4],
    rT[5], PET[1,5], PET[2,5], PET[3,5], PET[4,5], PET[5,5], PET[6,5],
    rT[6], PET[1,6], PET[2,6], PET[3,6], PET[4,6], PET[5,6], PET[6,6],
    rT[7], PET[1,7], PET[2,7], PET[3,7], PET[4,7], PET[5,7], PET[6,7],
    rT[8], PET[1,8], PET[2,8], PET[3,8], PET[4,8], PET[5,8], PET[6,8],
    rT[9], PET[1,9], PET[2,9], PET[3,9], PET[4,9], PET[5,9], PET[6,9],
    rT[10], PET[1,10], PET[2,10], PET[3,10], PET[4,10], PET[5,10], PET[6,10],
    rT[11], PET[1,11], PET[2,11], PET[3,11], PET[4,11], PET[5,11], PET[6,11],
    rT[12], PET[1,12], PET[2,12], PET[3,12], PET[4,12], PET[5,12], PET[6,12])
    
    return B_star_bar
    
end

function compute_Q_F1_SH1FC(aux)
    
    term1= divide_Vec12_into_Vec4Vec3(aux)
    term2 = compute_ssmatrices_Vec4Vec3_return_Vec4Mat33(term1)
    
    a = term2[1]
    b = term2[2]
    c = term2[3]
    d = term2[4]
    
    mat = Mat123(
    a[1,1], a[2,1], a[3,1], b[1,1], b[2,1], b[3,1], c[1,1], c[2,1], c[3,1], d[1,1], d[2,1], d[3,1], 
    a[1,2], a[2,2], a[3,2], b[1,2], b[2,2], b[3,2], c[1,2], c[2,2], c[3,2], d[1,2], d[2,2], d[3,2], 
    a[1,3], a[2,3], a[3,3], b[1,3], b[2,3], b[3,3], c[1,3], c[2,3], c[3,3], d[1,3], d[2,3], d[3,3])
    
    return mat
    
end 

function compute_D(D3)
    
    D311 = D3[1,1]
    D321 = D3[2,1]
    D331 = D3[3,1]
    
    D312 = D3[1,2]
    D322 = D3[2,2]
    D332 = D3[3,2]
    
    D313 = D3[1,3]
    D323 = D3[2,3]
    D333 = D3[3,3]
    
    E = Mat1212(
    D311, D321, D331, 0, 0, 0, -D311, -D321, -D331, 0, 0, 0,
    D312, D322, D332, 0, 0, 0, -D312, -D322, -D332, 0, 0, 0,
    D313, D323, D333, 0, 0, 0, -D313, -D323, -D333, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    -D311, -D321, -D331, 0, 0, 0, D311, D321, D331, 0, 0, 0,
    -D312, -D322, -D332, 0, 0, 0, D312, D322, D332, 0, 0, 0,
    -D313, -D323, -D333, 0, 0, 0, D313, D323, D333, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    
    return E
end

function compute_Bt(Tₛ⁻¹_th1g, Tₛ⁻¹_th2g)
    
    Bt = Mat1212(
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, Tₛ⁻¹_th1g[1,1], Tₛ⁻¹_th1g[2,1], Tₛ⁻¹_th1g[3,1], 0, 0, 0, 0, 0, 0,
    0, 0, 0, Tₛ⁻¹_th1g[1,2], Tₛ⁻¹_th1g[2,2], Tₛ⁻¹_th1g[3,2], 0, 0, 0, 0, 0, 0,
    0, 0, 0,  Tₛ⁻¹_th1g[1,3], Tₛ⁻¹_th1g[2,3], Tₛ⁻¹_th1g[3,3], 0, 0, 0, 0, 0, 0,
    0, 0, 0,  0, 0, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 0,  0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0,  0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, Tₛ⁻¹_th2g[1,1], Tₛ⁻¹_th2g[2,1], Tₛ⁻¹_th2g[3,1],
    0, 0, 0, 0, 0, 0, 0, 0, 0, Tₛ⁻¹_th2g[1,2], Tₛ⁻¹_th2g[2,2], Tₛ⁻¹_th2g[3,2],
    0, 0, 0, 0, 0, 0, 0, 0, 0, Tₛ⁻¹_th2g[1,3], Tₛ⁻¹_th2g[2,3], Tₛ⁻¹_th2g[3,3])
    
    return Bt
end

function compute_η_μ(Θ̄::AbstractVector{T}) where T

    Θ = norm(Θ̄)

    if Θ<10*eps(T)
        η = T(1/12)
        μ = T(1/360)
    else
        sinΘ, cosΘ = sincos(Θ)
        sinΘ2 = sin(Θ/2)
        η = ((2*sinΘ)-Θ*(1+cosΘ))/(2*(Θ^2)*sinΘ)
        μ = (Θ*(Θ+sinΘ)-8*(sinΘ2)^2)/(4*(Θ^4)*(sinΘ2)^2)
    end
    
    return η, μ

end



function compute_K̄ₕ(Θ̅, M̄, Tₛ⁻¹Θ̅)

    Θ̅M̄ᵀ = Θ̅*M̄'
    M̄Θ̅ᵀ = Θ̅M̄ᵀ'

    SΘ̅ = skew_skymmetric_matrix_from_vector(Θ̅)
    SM̄ = skew_skymmetric_matrix_from_vector(M̄)

    K̄ₕ =   η₁*(Θ̅M̄ᵀ - 2*M̄Θ̅ᵀ + dot(Θ̅, M̄)*ID3) * Tₛ⁻¹Θ̅
        + μ₁*(SΘ̅*SΘ̅*M̄*Θ̅') * Tₛ⁻¹Θ̅
        - (SM̄*Tₛ⁻¹Θ̅)/2

    return K̄ₕ

end





# rotation matrix from the global reference system to the rigidly_moving local reference
# system  of the deformed configuration
@inline function local_Rₑ(x₁, x₂, Rₑ⁰E₂)
    v₁ = (x₂ - x₁)/lₙ
    p₁ = R₁ * Rₑ⁰E₂
    p₂ = R₂ * Rₑ⁰E₂
    p = (p₁+p₂)/2
    v₃ = cross(v₁, p)
    v₃ = v₃ / norm(v₃)
    v₂ = cross(v₃, v₁)
    return [v₁ v₂ v₃]
end



@inline function auxiliary_variables(Rₑ, p, p₁, p₂)

    v₁ = Re[:,1]

    ru₁ = -v₁
    ru₂ = v₁

    q₁, q₂, _ = Rₑ*p
    q₁₁, q₁₂, _ = Rₑ*p₁
    q₂₁, q₂₂, _ = Rₑ*p₂
    
    η = q₁/q₂
    η₁₁ = q₁₁/q₂
    η₁₂ = q₁₂/q₂
    η₂₁ = q₂₁/q2
    η₂₂ = q₂₂/q2

    Gᵀu₁ = @Smatrix [0 0 η/lₙ; 0 0 1/lₙ; 0 -1/lₙ 0]
    GᵀΘ₁ = @Smatrix [η₁₂/2 -η₁₁/2 0; 0 0 0; 0 0 0]
    Gᵀu₂ = -Gᵀu₁
    GᵀΘ₂ = @Smatrix [η₂₂/2 -η₂₁/2 0; 0 0 0; 0 0 0]

    return ru₁, ru₂, η, (η₁₁, η₁₂, η₂₁, η₂₂), Gᵀu₁, GᵀΘ₁, Gᵀu₂, GᵀΘ₂

end



@inline function Pmatrices(N₁, N₂, N₃, N₄, N₅, N₆, lₙ, η, η₁₁, η₁₂, η₂₁, η₂₂)

    P₁Pu₁ = @Smatrix [0 0 0; 0 (N₃+N₄)/lₙ 0; 0 0 (N₃+N₄)/lₙ]
    P₁PΘ₁ = @Smatrix [0 0 0; 0 0 N₃; 0 -N₃ 0]
    P₁Pu₂ = -P₁Pu₁
    P₁PΘ₂ = @Smatrix [0 0 0; 0 0 N₄; 0 -N₄ 0]

    P₂Pu₁ = @Smatrix [0 0 -η*(N₁+N₂)/ln; 0 0 -(N₅+N₆)/lₙ; 0 (N₅+N₆)/lₙ 0]
    P₂PΘ₁ = @Smatrix [-(N₂*η₁₂)/2-N₁*(η₁₂/2 - 1) (η₁₁*(N₁+N₂))/2 0; 0 N₅ 0; 0 0 N₆]
    P₂Pu₂ = -P₂Pu₁
    P₂PΘ₂ = @Smatrix [-(N₁*η₂₂)/2-N₂*(η₂₂/2 - 1) (η₂₁*(N₁+N₂))/2 0; 0 N₅ 0; 0 0 N₆]


    return P₁Pu₁, P₁PΘ₁, P₁Pu₂, P₁PΘ₂, P₂Pu₁, P₂PΘ₁, P₂Pu₂, P₂PΘ₂

end




struct BlockMatDiag_4_4{T11, T22, T33, T44}
    B11::T11
    B22::T22
    B33::T33
    B44::T44
end

struct BlockMatDiag_2_2{T11, T22}
    B11::T11
    B22::T22
end

struct BlockMat_1_4{T11, T12, T13, T14}
    B11::T11
    B12::T12
    B13::T13
    B14::T14
end

struct BlockMat_2_4{T11, T12, T13, T14, T21, T22, T23, T24}
    B11::T11
    B12::T12
    B13::T13
    B14::T14
end

Base.:*(M1::BlockMat_1_4, M2::BlockMatDiag_4_4) = BlockMat_1_4(M1.B11*M2.B11, M1.B12*M2.B22, M1.B13*M2.B33, M1.B14*M2.B44)

Base.:*(M1::BlockMat_1_4, M2::BlockMat_4_4) = BlockMat_1_4(M1.B11*M2.B11 + M1.B12*M2.B21 + M1.B13*M2.B31 + M1.B14*M2.B41,
                                                           M1.B11*M2.B12 + M1.B12*M2.B22 + M1.B13*M2.B32 + M1.B14*M2.B42, 
                                                           M1.B11*M2.B13 + M1.B12*M2.B23 + M1.B13*M2.B33 + M1.B14*M2.B43, 
                                                           M1.B11*M2.B14 + M1.B12*M2.B24 + M1.B13*M2.B34 + M1.B14*M2.B44)