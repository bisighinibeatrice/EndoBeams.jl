test_sym(m) = maximum(abs.(Symmetric(m)-m))


# Compute beam contact contribution
function compute_contact_contribution(contact_vals, x1, x2, Re, comp, sdf, E, ddt, l0, EG, ln, fixed_matrices, rT, P, Tct, Kct, Cc, C_energy_e, iGP, e, sol_GP, to, T=Float64)
    
    @timeit_debug to "Get contact at the Gauss point" begin

        N1, N2, utG, H1G, wG, N7, P1G = contact_vals 
        xGv = N1*x1 + N2*x2 + Re*utG
        fc_eps, dfc_eps, Pic_eps, g_G, dg_G, ddg_G =  get_contact_GP(xGv, comp.epsC, sdf, T)

        sol_GP.xGP[e.indGP[iGP]] = xGv
        sol_GP.gGP[e.indGP[iGP]] = (g_G)/sdf.r
        sol_GP.status[e.indGP[iGP]] = 0
        sol_GP.fGP_N[e.indGP[iGP]] = zeros(Vec3{T})
        sol_GP.fGP_T[e.indGP[iGP]] = zeros(Vec3{T})

    end
    
    @timeit_debug to "Compute contact matrix and force" begin

        if fc_eps != 0

            # eq110 in [2]
            Inn = I - dg_G*dg_G'
            dgT = Inn*Re*H1G*E'*ddt
            dgTdgT = dot(dgT,dgT)
            
            G_eN = dg_G
            
            # eq112 in [2]
            G_eT = -comp.mu_T*dgT/sqrt(dgTdgT+comp.eps_tol_fric)
            
            sol_GP.fGP_N[e.indGP[iGP]] = fc_eps*G_eN
            sol_GP.fGP_T[e.indGP[iGP]] = fc_eps*G_eT
            sol_GP.status[e.indGP[iGP]] = 1

            # eqE2 in [2]
            G_e = G_eN + G_eT
            fC = fc_eps*G_e # fc_eps=pN(d) in [2]
            
            # ----------
            # fct_e term
            # ----------
            
            # eqE1 in [2]
            Tct += (l0/2)*wG*(E*H1G'*Re'*fC)
            
            # ----------
            # Kct term
            # ----------
            
            # term 1
            FC = Re'*fC  # eqE6 in [2] 
            aux = H1G'*FC
            
            SH1FC = compute_Q_F1_SH1FC(aux)
            
            term1 = -E*SH1FC*EG' # eqE18a in [2]
            
            # term 2
            SFC = get_skew_skymmetric_matrix_from_vector(FC)
            term2 = (N7/ln^2)*E*(fixed_matrices.A1)'*FC*rT - EG*SFC*P1G*P*E'
            
            # term 3
            term3 = E*H1G'*SFC*EG'
            
            # ----------
            #  term 4
            # ----------
            
            # term 4.1
            aux1 = E*H1G'*Re'*G_e
            aux2 = E*H1G'*Re'*dg_G
            term41 = dfc_eps*(aux1*aux2') # eqE16 in [2]
            
            # term 4.2
            
            # A1
            ReH1ET = Re*H1G*E'
            aux1 = dot(dg_G,ReH1ET*ddt)
            aux2 = ReH1ET'*ddg_G*ReH1ET*ddt
            term42_A1 = -aux1*(ddg_G*ReH1ET) - dg_G*aux2' # eqE18a in [2]
            
            # A2
            aux = H1G*E'*ddt
            SH1ETddt= get_skew_skymmetric_matrix_from_vector(aux)
            term42_A2 = -Inn*Re*SH1ETddt*EG' # eqE18b in [2]
            
            # A3
            aux = EG'*ddt
            SGEddt = get_skew_skymmetric_matrix_from_vector(aux)
            term42_A3 = Inn*(N7/(ln^2))*Re*fixed_matrices.A1*E'*ddt*rT + Inn*Re*SGEddt*P1G*P*E' # eqE18c in [2]
            
            # A4
            aux = E'*ddt         
            SDdot = compute_Q_F1_SH1FC(aux)
            
            term42_A4 = Inn*Re*H1G*SDdot*EG' # eqE14d in [2]
            
            term42_A = term42_A1 + term42_A2 + term42_A3 + term42_A4 # eqE13 in [2]
            
            aux1 = (I - (1/(dgTdgT+comp.eps_tol_fric))*(dgT*dgT'))*term42_A
            aux2 = ddg_G*Re*H1G*E'- comp.mu_T/(sqrt(dgTdgT+comp.eps_tol_fric))*aux1
            term42 = fc_eps*E*H1G'*Re'*aux2 # eqE16 in [2]
            
            term4 = term41 + term42  # eqE18d in [2]
            
            Kct += (l0/2)*wG*(term1 + term2 + term3 + term4)   # eqE19-22 in [2]
            
            # ----------
            # Cc term
            # ----------
            
            B = Inn*Re*H1G*E' # eqE14e in [2]
            aux1 = comp.mu_T/(sqrt(dgTdgT+comp.eps_tol_fric))
            aux = - fc_eps*aux1*E*H1G'*Re'*(I  - (1/(dgTdgT+comp.eps_tol_fric))*(dgT*dgT'))*B # eqE17 in [2]
            Cc += Cc + (l0/2)*wG*aux
            
            # energy
            C_energy_e = C_energy_e + (l0/2)*wG*Pic_eps
            
        end

    end
    
    return Tct, Kct, Cc, C_energy_e
    
end

# Compute beam dynamic contribution
function compute_dynamic_contribution(zG, l0, wG, P1G, P2G, NG, Θ1_bar, Θ2_bar, P, E, ddt, GT, ln, rT, Ddt_e, SWdt_e, Et, F1, ddtdt, Re, Tk, Tdamp, M, Ck, K_energy_e, fixed_matrices, conf, comp)
    
    xG = l0*(zG+1)/2
    
    # eqB3:5 in [2]: shape functions
    N1 = 1-xG/l0
    N2 = 1-N1
    N3 = xG*(1-xG/l0)^2
    N4 = -(1-xG/l0)*((xG^2)/l0)
    N5 = (1-3*xG/l0)*(1-xG/l0)
    N6 = (3*xG/l0-2)*(xG/l0)
    N7 = N3+N4
    N8 = N5+N6-1
    
    # eqD4 in [2]: local transerve displacement vector
    utG = reshapeP1G(N3, N4, Θ1_bar, Θ2_bar)
    SutG = get_skew_skymmetric_matrix_from_vector(utG)
    
    # eqD14 in [2]
    utdtG = P1G*P * E' * ddt
    SutdtG = get_skew_skymmetric_matrix_from_vector(utdtG)
    
    # ???: local transerve rotational vector
    ΘG_bar = P2G * Vec6(Θ1_bar[1], Θ1_bar[2], Θ1_bar[3],  Θ2_bar[1], Θ2_bar[2], Θ2_bar[3])
    SuΘG_bar = get_skew_skymmetric_matrix_from_vector(ΘG_bar)
    
    # eq74 in [2]: RG_bar, rotation matrix associated with the local cross-section rotation
    RG_bar = ID3 + SuΘG_bar
    
    # eq100 in [2]
    IrhoeG = RG_bar*(conf.mat.Jrho)*RG_bar'
    
    # ------------------------------------
    # expression of H1, H2 and H1dt, H2dt
    
    # eqD6 in [2], matrix H1
    H1G = NG + P1G*P - SutG*GT
    
    # eqD16 in [2], matrix H2
    H2G = P2G*P + GT
    
    # eqD14 in [2], matrix H1dt
    H1dtG = (N7/(ln^2))*fixed_matrices.A1*(rT*ddt)-SutdtG*GT
    
    # eqD17 in [2], matrix H2dt
    H2dtG = (N8/(ln^2))*fixed_matrices.A2*(rT*ddt)
    
    # -----------------------------
    # expression of C1, C2, C3, C4
    
    # eqD33 in [2], h1G
    h1G = Vec3(H1G * Ddt_e)
    
    # eqD34 in [2], h2G
    h2G = Vec3(H2G * Ddt_e)
    
    # skew-symmetric matrix of h1
    Sh1G = get_skew_skymmetric_matrix_from_vector(h1G)
    
    # skew-symmetric matrix of h2
    Sh2G = get_skew_skymmetric_matrix_from_vector(h2G)
    
    # eqD13 in [2], C1
    C1G = SWdt_e*H1G + H1dtG - H1G*Et
    
    # eqD22 in [2], C2
    C2G = SWdt_e*H2G + H2dtG - H2G*Et
    
    # eqD30 in [2], C3
    C3G = -Sh1G*GT + N7/(ln^2)*fixed_matrices.A1*Ddt_e*(fixed_matrices.re)' + SWdt_e*P1G*P + H1G*F1*GT

    # eqD31 in [2], C4
    C4G = -Sh2G*GT + N8/(ln^2)*fixed_matrices.A2*Ddt_e*(fixed_matrices.re)' + H2G*F1*GT
    
    # --------------------------------
    # expression of Udtdt, Wdt, Wdtdt
    
    # eqD12 in [2], udtdt --> Re
    Udtdt = H1G*E'*ddtdt + C1G*E'*ddt
    
    # eqD20 in [2], wdtdt --> Re'
    Wdt = H2G*E'*ddt
    SWdt = get_skew_skymmetric_matrix_from_vector(Wdt)
    
    # eqD21 in [2], wdtdt --> Re'
    Wdtdt = H2G*E'*ddtdt + C2G*E'*ddt
    
    # ---------------------------------
    # expression of Tk_e, M_e and Ck_e
    
    aux1 = (l0/2)*wG*conf.mat.Arho*H1G'
    wGH2GT = (l0/2)*wG*H2G'
    aux2 = wGH2GT*IrhoeG
    aux3 = wGH2GT*SWdt
    
    Udt = H1G*E'*ddt
    
    term1_Tk = aux1*Udtdt
    term2_Tk = aux2*Wdtdt
    term3_Tk = aux3*IrhoeG*Wdt
    term_viscosity_1_Tk = comp.damping*aux1*Udt
    term_viscosity_2_Tk = comp.damping*aux2*Wdt
    
    # eq98 in [2]
    Tk += term1_Tk + term2_Tk + term3_Tk + term_viscosity_1_Tk + term_viscosity_2_Tk
    Tdamp += term_viscosity_1_Tk + term_viscosity_2_Tk
    
    # eq102 in [2]
    M +=  aux1*H1G + aux2*H2G

    
    term1_Ck = aux1 * (C1G + C3G)
    term2_Ck = aux2 * (C2G + C4G)
    IrhoeGWdt = IrhoeG*Wdt
    SIrhoeGWdt = get_skew_skymmetric_matrix_from_vector(IrhoeGWdt)
    term3_Ck = (aux3*IrhoeG)*H2G
    term4_Ck = -(wGH2GT*SIrhoeGWdt)*H2G
    term_viscosity_1_Ck = comp.damping*aux1*H1G
    term_viscosity_2_Ck = comp.damping*aux2*H2G
    
    # eq103 in [2]
    Ck += term1_Ck + term2_Ck + term3_Ck + term4_Ck + term_viscosity_1_Ck + term_viscosity_2_Ck
    
    # kinetic energy
    IrhoG = Re*IrhoeG*Re'
    wdt_G = Re*Wdt
    udt_eG =  H1G*E'*ddt
    udt_G = Re*udt_eG
    
    K_energy_e = K_energy_e + (l0/2)*wG/2*(conf.mat.Arho*(udt_G'*udt_G)) + (l0/2)*wG/2*(wdt_G'*(IrhoG*wdt_G))
    
    return Tk, Tdamp, M, Ck, K_energy_e, (N1, N2, utG, H1G, wG, N7, P1G)
    
end 

# Compute beam internal energy
function compute_internal_energy(conf, D_italic_bar, l0)
    
    mat = conf.mat
    geom = conf.geom
    
    A = geom.A
    Emat = mat.E
    Gmat = mat.G
    I22 = geom.I22
    I33 = geom.I33
    J = geom.J
    
    term1_energy = A*Emat*(D_italic_bar[1]^2)/2
    term2_energy = Emat*I22*(D_italic_bar[4]^2)/2
    term3_energy = Emat*I22*(D_italic_bar[7]^2)/2
    term4_energy = Emat*I33*(D_italic_bar[3]^2)/2
    term5_energy = Emat*I33*(D_italic_bar[6]^2)/2
    term6_energy = Gmat*J*(D_italic_bar[2]^2)/2
    term7_energy = Gmat*J*(D_italic_bar[5]^2)/2
    term8_energy = Emat*I33*D_italic_bar[3]*D_italic_bar[6]/2
    term9_energy = Emat*I22*D_italic_bar[4]*D_italic_bar[7]/2
    term10_energy = - Gmat*J*D_italic_bar[2]*D_italic_bar[5]
    
    Phi_energy_e = (term1_energy+term2_energy+term3_energy+term4_energy+term5_energy+term6_energy+term7_energy+term8_energy+term9_energy+term10_energy)/l0
    
    return Phi_energy_e
    
end

# Compute beam contributions
function get_beam_contributions!(e, allnodes, conf, sdf, fixed_matrices, comp, sol_GP, to, T=Float64) 
    
    @timeit_debug to "Move to local reference system" begin
        
        # retrieve the matrix Re_0 of the beam
        Re0 = e.R0
        
        # information from node 1 and 2
        X1, X2 = get_local_pos(e, allnodes)
        uk1, uk2 = get_local_displ(e, allnodes)
        ukdt1, ukdt2 = get_local_vel(e, allnodes)
        ukdtdt1, ukdtdt2 = get_local_acc(e, allnodes)
        wkdt1, wkdt2 = get_local_ang_vel(e, allnodes)
        wkdtdt1, wkdtdt2 = get_local_ang_acc(e, allnodes)
        R1, R2 = get_local_rot(e, allnodes)
        
        # current position of node i1 (global rs)
        x1 =  X1 + uk1
        x2 =  X2 + uk2
        
        # -------------------------------------------------------------------------------------------
        # rigidly_moving rs of the deformed configuration (v1, v2, v3)
        # -------------------------------------------------------------------------------------------
        # with respect to the new rs, the first node is placed at the origin (0) and the second at a
        # distance l: i1 --> 0 and i2 --> l
        
        # length of the beam in the undeformed configuration (CHECK)
        l0 = e.l0
        
        #  length of the beam in the deformed configuration
        ln = norm(x2 - x1)
        
        # eq32a in [2]: local displacement with respect to the local rs of the deformed configuration
        # (ux_0=uy_0=uz_0=uyL=uz_L=0 as a consequence of how the local rs (v1, v2, v3) is built)
        u_bar = ln - l0
        
        #  eq23 in [2]: v1
        v1 = (x2 - x1)/ln
        
        #  eq24 in [2]: auxiliary vector p
        t2_0 = R1 * Re0 * Vec3(0,1,0)
        t2_l = R2 * Re0 * Vec3(0,1,0)
        p = (t2_0 + t2_l)/2
        
        #  eq23 in [2]: v2 and v3
        v3 = cross(v1, p)
        v3 = v3 / norm(v3)
        v2 = cross(v3, v1)
        
        # rotation matrix from the global reference system to the rigidly_moving local reference
        # system  of the deformed configuration
        Re =  [v1 v2 v3]
        
        # eq70 in [2]: matrix E = diag(Re)
        E = compute_E_or_Et(Re)
        
        # eqC4 in [2]: r
        r = compute_r(v1)
        rT = r' #TODO
        
        # -------------------------------------------------------------------------------------------
        # D_italic_bar = [u_bar, Θ1_bar, Θ2_bar]
        # -------------------------------------------------------------------------------------------
        
        # eq72:73 in [2]: (local rotations) rotation matrix from the deformation-
        # undependent local rs to the deformation-dependent local rs of the
        # deformed configuration
        R1_bar = Re' * R1 * Re0
        R2_bar = Re' * R2 * Re0
        
        # Theta1_bar
        Θ1_bar = get_angle_from_rotation_matrix(R1_bar)
        
        # Theta2_bar
        Θ2_bar = get_angle_from_rotation_matrix(R2_bar)
        
        # eq78 in [2]: D_italic_bar
        D_italic_bar = Vec7(u_bar, Θ1_bar[1], Θ1_bar[2], Θ1_bar[3], Θ2_bar[1], Θ2_bar[2], Θ2_bar[3])
        
    end 
    
    @timeit_debug to "Compute internal forces vector, stiffness matrix" begin
        
        # -------------------------------------------------------------------------------------------
        # -------------------------------------------------------------------------------------------
        # Compute internal forces vector, stiffness matrix
        # -------------------------------------------------------------------------------------------
        # -------------------------------------------------------------------------------------------
        
        # -------------------------------------------------------------------------------------------
        # Kint_bar and Fint_bar
        Kint_bar = e.Kint
        Fint_bar = Kint_bar * D_italic_bar
        
        # -------------------------------------------------------------------------------------------
        # Internal energy    
        
        Phi_energy_e = compute_internal_energy(conf, D_italic_bar, l0)
        
        # -------------------------------------------------------------------------------------------
        # from Kint_bar and Fint_bar to Kint and Fint
        # -------------------------------------------------------------------------------------------
        
        # inverse of the screw_symmetric matrix of Θ1_bar and Θ2_bar
        TsinvΘ1_bar = get_inverse_skew_skymmetric_matrix_from_angle(Θ1_bar)
        TsinvΘ2_bar = get_inverse_skew_skymmetric_matrix_from_angle(Θ2_bar)
        
        # eqC2 in [2]: B_bar
        B_bar = compute_B_bar(TsinvΘ1_bar, TsinvΘ2_bar) 
        
        GT, P, eta = compute_GT_P_eta(ln, Re, p, t2_0,  t2_l)
        G = GT'
        EG = E*G
        
        #  eqC3 in [2]: B_star_bar
        B_star_bar = compute_B_star_bar(rT, P*E') 
        
        # eqC8 in [2]: B
        B = B_bar * B_star_bar 
        
        # ------------------
        # eq90 in [2]: Fint
        # ------------------
        
        Tint = B' * Fint_bar
        
        # eqC12/13 in [2]:  Kha_bar
        M1 = Vec3(Fint_bar[2], Fint_bar[3], Fint_bar[4])
        M2 = Vec3(Fint_bar[5], Fint_bar[6], Fint_bar[7])
        Kh1_bar = get_Kah_bar(Θ1_bar, M1, TsinvΘ1_bar)
        Kh2_bar = get_Kah_bar(Θ2_bar, M2, TsinvΘ2_bar)
        
        # eqC11 in [2]:  Kh_bar
        Kh_bar = compute_Kh_bar(Kh1_bar, Kh2_bar) 
        
        # eqC18 in [2]: F_star_bar = [N_bar m]
        F_star_bar = B_bar' * Fint_bar
        
        # eqC17 in [2]: PTm = [Q1 Q2 Q3 Q4]
        m = Vec6(F_star_bar[2], F_star_bar[3], F_star_bar[4], F_star_bar[5], F_star_bar[6], F_star_bar[7])
        PTm = P' * m
        
        # eqC16 in [2]: Q
        Q = compute_Q_F1_SH1FC(PTm)
        
        # eqC16 in [2]: a
        a = Vec3(0, eta*(m[1]+m[4])/ln - (m[2]+m[5])/ln, (m[3]+m[6])/ln)
        
        # eqC15 in [2]: D3
        D3 = (ID3 - v1*v1')/ln
        
        # eqC15 in [2]: D
        D = compute_D(D3)
        
        # eqC14 in [2]: Km
        Km = D * F_star_bar[1] - E*Q*GT*E' + EG*a*rT 
        
        # -------------------
        # eqC10 in [2]: Kint
        # -------------------
        
        Kint = B_star_bar'*( B_bar'*Kint_bar*B_bar + Kh_bar)*B_star_bar + Km 
    end

    @timeit_debug to "Dynamic and contact contributions" begin

        
        # Gauss quadrature
        nG = comp.nG
        wG_v = comp.wG
        zG_v = comp.zG
        
        # eq69-D32 in [2]: d = [u1, θ1, u2, θ2], ddt = [udt1, wdt1, udt2, wdt2], ...
        ddt = Vec12(ukdt1[1], ukdt1[2], ukdt1[3], wkdt1[1], wkdt1[2], wkdt1[3], ukdt2[1], ukdt2[2], ukdt2[3], wkdt2[1], wkdt2[2], wkdt2[3])
        ddtdt = Vec12(ukdtdt1[1], ukdtdt1[2], ukdtdt1[3], wkdtdt1[1], wkdtdt1[2], wkdtdt1[3], ukdtdt2[1], ukdtdt2[2], ukdtdt2[3], wkdtdt2[1], wkdtdt2[2], wkdtdt2[3])
        
        # eq70 in [2]: D_italic = [U1, Θ1, U2, Θ2], Ddt_italic = [Udt1, Wdt1, Udt2, Wdt2]
        Ddt_e = E' * ddt
        
        # skew-symmetric matrices of velocity Udt1 and Udt1 angular velocity
        F1 = compute_Q_F1_SH1FC(Ddt_e)
        
        # eqD3 in [2]: Wdt ???
        Wdt_e = GT * E' * ddt
        
        # skew-symmetric matrices of Wdt ???
        SWdt_e = get_skew_skymmetric_matrix_from_vector(Wdt_e)
        
        # reference?? Et=E appendix D in [2]
        Et = compute_E_or_Et(SWdt_e)
        
        # initialise the local matrices used in the Gauss loop
        K_energy_e = zero(T)
        C_energy_e = zero(T)
        
        Tk = zeros(Vec12{T})
        Tdamp = zeros(Vec12{T})
        M = zeros(Mat1212{T}) 
        Ck = zeros(Mat1212{T})
        Tct = zeros(Vec12{T}) 
        Kct = zeros(Mat1212{T}) 
        Cc = zeros(Mat1212{T}) 

        # cycle among the Gauss positions
        for iG in 1:nG

            @timeit_debug to "Dynamic contributions" begin

                Tk, Tdamp, M, Ck, K_energy_e, contact_vals = compute_dynamic_contribution(zG_v[iG], l0,  wG_v[iG], fixed_matrices.P1G_v[iG], fixed_matrices.P2G_v[iG], fixed_matrices.NG_v[iG], Θ1_bar, Θ2_bar, P, E, ddt, GT, ln, rT, Ddt_e, SWdt_e, Et, F1, ddtdt, Re, Tk, Tdamp, M, Ck, K_energy_e, fixed_matrices, conf, comp)
            end 

            @timeit_debug to "Contact contributions" begin

                if !isnothing(sdf)
                    Tct, Kct, Cc, C_energy_e =  compute_contact_contribution(contact_vals, x1, x2, Re, comp, sdf, E, ddt, l0, EG, ln, fixed_matrices, rT, P, Tct, Kct, Cc, C_energy_e, iG, e, sol_GP, to, T)
                end     

            end 
            
        end 

    end 
         
    @timeit_debug to "Return to global reference system" begin
        
        Tk = E * Tk
        Tdamp =  E * Tdamp
        M = E * M * E'
        Ck = E * Ck * E'
        

        # Note: this has no effect on convergence rate from my tests
        
        # Deltn1_1, Deltn1_2 = get_local_rot_delt(e, allnodes)
        
        # thi_n1_1 = get_angle_from_rotation_matrix(Deltn1_1)
        # thi_n1_2 = get_angle_from_rotation_matrix(Deltn1_2)
        
        # Tsinv_th1g = get_inverse_skew_skymmetric_matrix_from_angle(thi_n1_1)
        # Tsinv_th2g = get_inverse_skew_skymmetric_matrix_from_angle(thi_n1_2)
        
        # Bt = compute_Bt(Tsinv_th1g, Tsinv_th2g)
        
        # M = M*Bt
        # Ck = Ck*Bt
        # Cc = Cc*Bt
        
        Ck += Cc 
        
    end 

    
    return (Kint, Tint, Tk, Tdamp, Ck, M, Tct, Kct), (Phi_energy_e, K_energy_e, C_energy_e)
    
end 

# Compute elemental contributions and assembly
function compute_K_T!(allnodes, allbeams, matrices, energy, conf, sdf, fixed_matrices, comp, sol_GP, to, T=Float64) 
    
    @timeit_debug to "Initialization" begin
        
        # initialise the matrices associate to the whole structure
        matrices.Kint.nzval .= 0
        matrices.Ck.nzval .= 0 
        matrices.M.nzval .= 0
        matrices.Kct.nzval .= 0
        matrices.Tint .= 0
        matrices.Tk .= 0
        matrices.Tdamp .= 0
        matrices.Tct .= 0
       
        # initialise the energy values associate to the whole structure
        energy.Phi_energy = 0
        energy.K_energy = 0
        energy.C_energy = 0
        
    end
    
    # cycle over the beams   
    @timeit_debug to "Compute and assemble elemental contributions" begin

        for e in LazyRows(allbeams)
            
            @timeit_debug to "Compute elemental contributions" begin

                #----------------------------------------
                # Compute the contibution from the e beam
                (Kint, Tint, Tk, Tdamp, Ck, M, Tct, Kct), (Phi_energy, K_energy, C_energy)  = get_beam_contributions!(e, allnodes, conf, sdf, fixed_matrices, comp, sol_GP, to, T) 
            
            end 

            @timeit_debug to "Assemble elemental contributions" begin

                #-----------------------
                # Assemble contributions
                
                i1 = e.node1
                i2 = e.node2
                idof1 = allnodes.idof_6[i1]
                idof2 = allnodes.idof_6[i2]
                
                idof = vcat(idof1, idof2)
                
                energy.Phi_energy +=  Phi_energy
                energy.K_energy += K_energy
                energy.C_energy +=  C_energy
            
                update_spmat_sum(matrices.Kint, e.sparsity_map, Kint)
                update_spmat_sum(matrices.Ck, e.sparsity_map, Ck)
                update_spmat_sum(matrices.M, e.sparsity_map, M)
                update_spmat_sum(matrices.Kct, e.sparsity_map, Kct)

                
                update_vec_sum(matrices.Tk, idof, Tk)
                update_vec_sum(matrices.Tdamp, idof, Tdamp)
                update_vec_sum(matrices.Tint, idof, Tint)
                update_vec_sum(matrices.Tct, idof, Tct)
            end
                               
        end

        println("Kint: $(test_sym(matrices.Kint))")
        println("Ck: $(test_sym(matrices.Ck))")
        println("M: $(test_sym(matrices.M))")
        println("Kct: $(test_sym(matrices.Kct))")



    end 


    
end 