#-------------------------------------------------------------------
# FUNCTIONS AT THE BEGINNING OF EACH BEAM CONTRIBUTIONS COMPUTATION
#-------------------------------------------------------------------

# Gets the position of the two nodes forming the given beam element
function get_local_pos(e, allnodes)
    
    return allnodes.pos[e.node1], allnodes.pos[e.node2]
    
end 

# Gets the displacement of the two nodes forming the given beam element
function get_local_displ(e, allnodes)
    
    return allnodes.u[e.node1], allnodes.u[e.node2]
    
end 

# Gets the velocity of the two nodes forming the given beam element
function get_local_vel(e, allnodes)
    
    return allnodes.udt[e.node1], allnodes.udt[e.node2]
    
end 

# Gets the acceleration of the two nodes forming the given beam element
function get_local_acc(e, allnodes)
    
    return allnodes.udtdt[e.node1], allnodes.udtdt[e.node2]
    
end 

# Gets the angular velocity of the two nodes forming the given beam element
function get_local_ang_vel(e, allnodes)
    
    return allnodes.wdt[e.node1], allnodes.wdt[e.node2]
    
end 

# Gets the angular acceleration of the two nodes forming the given beam element
function get_local_ang_acc(e, allnodes)
    
    return allnodes.wdtdt[e.node1], allnodes.wdtdt[e.node2]
    
end 

# Gets the rotation matrix of the two nodes forming the given beam element
function get_local_rot(e, allnodes)
    
    return allnodes.R[e.node1], allnodes.R[e.node2]
    
end 

# Gets the rotation matrix displacement of the two nodes forming the given beam element
function get_local_rot_delt(e, allnodes)
    
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
    
    aux1 = get_skew_skymmetric_matrix_from_vector(a[1])
    aux2 = get_skew_skymmetric_matrix_from_vector(a[2])
    aux3 = get_skew_skymmetric_matrix_from_vector(a[3])
    aux4 = get_skew_skymmetric_matrix_from_vector(a[4])
    
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

function reshapeP1G(N3, N4, Θ1_bar, Θ2_bar)
    
    return Vec3(0, N3*Θ1_bar[3]+N4*Θ2_bar[3], -(N3*Θ1_bar[2]+N4*Θ2_bar[2]))
    
end 

#------------------------------------------
# FUNCTIONS TO COMPUTE SPECIFIC MATRICES
#------------------------------------------

function compute_E_or_Et(Re)
    
    Re11 = Re[1,1]
    Re21 = Re[2,1]
    Re31 = Re[3,1]
    
    Re12 = Re[1,2]
    Re22 = Re[2,2]
    Re32 = Re[3,2]
    
    Re13 = Re[1,3]
    Re23 = Re[2,3]
    Re33 = Re[3,3]
    
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
    
    return Symmetric(E)
end

function compute_r(v)
    
    r = Vec12(-v[1], -v[2], -v[3], 0, 0, 0, v[1], v[2], v[3], 0, 0, 0)
    
    return r
    
end

function compute_GT_P_eta(ln, Re, p, t2_0,  t2_l)
    
    # eqC7 in [2]: eta, eta11, eta12, eta21, eta22
    aux = Re'*p; q1 = aux[1]; q2 = aux[2]
    aux = Re'*t2_0; q11 = aux[1]; q12 = aux[2]
    aux = Re'*t2_l; q21 = aux[1]; q22 = aux[2]
    
    eta = q1/q2
    eta11 = q11/q2
    eta12 = q12/q2
    eta21 = q21/q2
    eta22 = q22/q2
    
    # eqC6 in [2]: GT
    GT = Mat312(
    0, 0, 0,
    0, 0, -1/ln,
    eta/ln, 1/ln, 0,
    eta12/2, 0, 0,
    -eta11/2, 0, 0,
    0, 0, 0,
    0, 0, 0,
    0, 0, 1/ln,
    -eta/ln, -1/ln, 0,
    eta22/2, 0, 0,
    -eta21/2, 0, 0,
    0, 0, 0)
    
    P = Mat612(
    0, 0, 0, 0, 0, 0,
    0, 0, 1/ln, 0, 0, 1/ln,
    -eta/ln, -1/ln, 0, -eta/ln, -1/ln, 0,
    1-eta12/2, 0, 0, -eta12/2, 0, 0,
    eta11/2, 1, 0, eta11/2, 0, 0,
    0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0,
    0, 0, -1/ln, 0, 0, -1/ln,
    eta/ln, 1/ln, 0, eta/ln, 1/ln, 0,
    -eta22/2, 0, 0, 1-eta22/2, 0, 0,
    eta21/2, 0, 0, eta21/2, 1, 0,
    0, 0, 0, 0, 0, 1)
    
    return GT, P, eta
    
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

function compute_B_bar(TsinvΘ1_bar, TsinvΘ2_bar) 
    
    B_bar = Mat77(
    1, 0, 0, 0, 0, 0, 0,
    0, TsinvΘ1_bar[1,1], TsinvΘ1_bar[2,1], TsinvΘ1_bar[3,1], 0, 0, 0,
    0, TsinvΘ1_bar[1,2], TsinvΘ1_bar[2,2], TsinvΘ1_bar[3,2], 0, 0, 0,
    0, TsinvΘ1_bar[1,3], TsinvΘ1_bar[2,3], TsinvΘ1_bar[3,3], 0, 0, 0,
    0, 0, 0, 0, TsinvΘ2_bar[1,1], TsinvΘ2_bar[2,1], TsinvΘ2_bar[3,1],
    0, 0, 0, 0, TsinvΘ2_bar[1,2], TsinvΘ2_bar[2,2], TsinvΘ2_bar[3,2],
    0, 0, 0, 0, TsinvΘ2_bar[1,3], TsinvΘ2_bar[2,3], TsinvΘ2_bar[3,3])
    
    return B_bar
    
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

function compute_Bt(Tsinv_th1g, Tsinv_th2g)
    
    Bt = Mat1212(
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, Tsinv_th1g[1,1], Tsinv_th1g[2,1], Tsinv_th1g[3,1], 0, 0, 0, 0, 0, 0,
    0, 0, 0, Tsinv_th1g[1,2], Tsinv_th1g[2,2], Tsinv_th1g[3,2], 0, 0, 0, 0, 0, 0,
    0, 0, 0,  Tsinv_th1g[1,3], Tsinv_th1g[2,3], Tsinv_th1g[3,3], 0, 0, 0, 0, 0, 0,
    0, 0, 0,  0, 0, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 0,  0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0,  0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, Tsinv_th2g[1,1], Tsinv_th2g[2,1], Tsinv_th2g[3,1],
    0, 0, 0, 0, 0, 0, 0, 0, 0, Tsinv_th2g[1,2], Tsinv_th2g[2,2], Tsinv_th2g[3,2],
    0, 0, 0, 0, 0, 0, 0, 0, 0, Tsinv_th2g[1,3], Tsinv_th2g[2,3], Tsinv_th2g[3,3])
    
    return Bt
end

function get_eta_mu(theta::AbstractVector{T}) where T
    
    eta_1 = 1/12
    mu_1 = 1/360
    
    theta= norm(theta)
    
    sin_theta = sin(theta)
    cos_theta = cos(theta)
    sin_theta2 = sin(theta/2)
    
    eta_2 = ((2*sin_theta)-theta*(1+cos_theta))/(2*(theta^2)*sin_theta)
    
    mu_2 = (theta*(theta+sin_theta)-8*(sin_theta2)^2)/(4*(theta^4)*(sin_theta2)^2)
    
    ind_choose = theta<10*eps(T) ? 1 : 2
    
    aux_eta = (eta_1, eta_2)
    aux_mu = (mu_1, mu_2)
    
    
    return aux_eta[ind_choose], aux_mu[ind_choose]
    
end

function get_Kah_bar(theta, M, Tsinvtheta)
    
    Stheta = get_skew_skymmetric_matrix_from_vector(theta)
    SM = get_skew_skymmetric_matrix_from_vector(M)
    
    eta, mu = get_eta_mu(theta)
    
    aux = theta*M'
    
    aux1 = eta * (aux - 2*aux'+ (theta'*M)*ID3) 
    aux2 = mu * (Stheta * Stheta * aux') 
    aux3 = - SM/2
    
    Kah_bar = (aux1 + aux2 + aux3) * Tsinvtheta
    
    return Kah_bar
    
end 


