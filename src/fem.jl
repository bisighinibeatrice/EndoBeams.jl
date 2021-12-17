
# Compute beam contact contribution
@fastmath function compute_contact_contribution(contact_vals, x₁, x₂, Rₑ, comp, sdf, E, ddt, l₀, EG, lₙ, fixed_matrices, rT, P, Tct, Kct, Cc, C_energy_e, iGP, e, sol_GP, T=Float64)
    
    # @timeit_debug to "Get contact at the Gauss point" begin

        N1, N2, utG, H1G, ωG, N7, P1G = contact_vals 
        xGv = N1*x₁ + N2*x₂ + Rₑ*utG
        fc_eps, dfc_eps, Pic_eps, g_G, dg_G, ddg_G =  get_contact_GP(xGv, comp.epsC, sdf, T)

        sol_GP.xGP[e.indGP[iGP]] = xGv
        sol_GP.gGP[e.indGP[iGP]] = (g_G)/sdf.r
        sol_GP.status[e.indGP[iGP]] = 0
        sol_GP.fGP_N[e.indGP[iGP]] = zeros(Vec3{T})
        sol_GP.fGP_T[e.indGP[iGP]] = zeros(Vec3{T})

    # end
    
    # @timeit_debug to "Compute contact matrix and force" begin

        if fc_eps != 0

            # eq110 in [2]
            Inn = ID3 - dg_G*dg_G'
            dgT = Inn*Rₑ*H1G*E'*ddt
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
            Tct += (l₀/2)*ωG*(E*H1G'*Rₑ'*fC)
            
            # ----------
            # Kct term
            # ----------
            
            # term 1
            FC = Rₑ'*fC  # eqE6 in [2] 
            aux = H1G'*FC
            
            SH1FC = compute_Q_F1_SH1FC(aux)
            
            term1 = -E*SH1FC*EG' # eqE18a in [2]
            
            # term 2
            SFC = skew_skymmetric_matrix_from_vector(FC)
            term2 = (N7/lₙ^2)*E*(fixed_matrices.A1)'*FC*rT - EG*SFC*P1G*P*E'
            
            # term 3
            term3 = E*H1G'*SFC*EG'
            
            # ----------
            #  term 4
            # ----------
            
            # term 4.1
            aux1 = E*H1G'*Rₑ'*G_e
            aux2 = E*H1G'*Rₑ'*dg_G
            term41 = dfc_eps*(aux1*aux2') # eqE16 in [2]
            
            # term 4.2
            
            # A1
            ReH1ET = Rₑ*H1G*E'
            aux1 = dot(dg_G,ReH1ET*ddt)
            aux2 = ReH1ET'*ddg_G*ReH1ET*ddt
            term42_A1 = -aux1*(ddg_G*ReH1ET) - dg_G*aux2' # eqE18a in [2]
            
            # A2
            aux = H1G*E'*ddt
            SH1ETddt= skew_skymmetric_matrix_from_vector(aux)
            term42_A2 = -Inn*Rₑ*SH1ETddt*EG' # eqE18b in [2]
            
            # A3
            aux = EG'*ddt
            SGEddt = skew_skymmetric_matrix_from_vector(aux)
            term42_A3 = Inn*(N7/(lₙ^2))*Rₑ*fixed_matrices.A1*E'*ddt*rT + Inn*Rₑ*SGEddt*P1G*P*E' # eqE18c in [2]
            
            # A4
            aux = E'*ddt         
            SDdot = compute_Q_F1_SH1FC(aux)
            
            term42_A4 = Inn*Rₑ*H1G*SDdot*EG' # eqE14d in [2]
            
            term42_A = term42_A1 + term42_A2 + term42_A3 + term42_A4 # eqE13 in [2]
            
            aux1 = (ID3 - (1/(dgTdgT+comp.eps_tol_fric))*(dgT*dgT'))*term42_A
            aux2 = ddg_G*Rₑ*H1G*E'- comp.mu_T/(sqrt(dgTdgT+comp.eps_tol_fric))*aux1
            term42 = fc_eps*E*H1G'*Rₑ'*aux2 # eqE16 in [2]
            
            term4 = term41 + term42  # eqE18d in [2]
            
            Kct += (l₀/2)*ωG*(term1 + term2 + term3 + term4)   # eqE19-22 in [2]
            
            # ----------
            # Cc term
            # ----------
            
            B = Inn*Rₑ*H1G*E' # eqE14e in [2]
            aux1 = comp.mu_T/(sqrt(dgTdgT+comp.eps_tol_fric))
            aux = - fc_eps*aux1*E*H1G'*Rₑ'*(ID3  - (1/(dgTdgT+comp.eps_tol_fric))*(dgT*dgT'))*B # eqE17 in [2]
            Cc += Cc + (l₀/2)*ωG*aux
            
            # energy
            C_energy_e = C_energy_e + (l₀/2)*ωG*Pic_eps
            
        end

    # end
    
    return Tct, Kct, Cc, C_energy_e
    
end

# Compute beam dynamic contribution
@fastmath function compute_dynamic_contribution(zG, l₀, ωG, P1G, P2G, NG, Θ̅₁, Θ̅₂, P, E, ddt, GT, lₙ, rT, Ddt_e, SWdt_e, Et, F1, ddtdt, Rₑ, Tk, Tdamp, M, Ck, K_energy_e, fixed_matrices, conf, comp)
    
    ξ = l₀*(zG+1)/2
    
    # eqB3:5 in [2]: shape functions
    N1 = 1-ξ/l₀
    N2 = 1-N1
    N3 = ξ*(1-ξ/l₀)^2
    N4 = -(1-ξ/l₀)*((ξ^2)/l₀)
    N5 = (1-3*ξ/l₀)*(1-ξ/l₀)
    N6 = (3*ξ/l₀-2)*(ξ/l₀)
    N7 = N3+N4
    N8 = N5+N6-1
    
    # eqD4 in [2]: local transerve displacement vector
    utG = reshapeP1G(N3, N4, Θ̅₁, Θ̅₂)
    SutG = skew_skymmetric_matrix_from_vector(utG)
    
    # eqD14 in [2]
    utdtG = P1G*P * E' * ddt
    SutdtG = skew_skymmetric_matrix_from_vector(utdtG)
    
    # ???: local transerve rotational vector
    ΘG_bar = P2G * Vec6(Θ̅₁[1], Θ̅₁[2], Θ̅₁[3],  Θ̅₂[1], Θ̅₂[2], Θ̅₂[3])
    SuΘG_bar = skew_skymmetric_matrix_from_vector(ΘG_bar)
    
    # eq74 in [2]: RG_bar, rotation matrix associated with the local cross-section rotation
    RG_bar = ID3 + SuΘG_bar
    
    # eq100 in [2]
    IrhoeG = RG_bar*(conf.mat.Jᵨ)*RG_bar'
    
    # ------------------------------------
    # expression of H1, H2 and H1dt, H2dt
    
    # eqD6 in [2], matrix H1
    H1G = NG + P1G*P - SutG*GT
    
    # eqD16 in [2], matrix H2
    H2G = P2G*P + GT
    
    # eqD14 in [2], matrix H1dt
    H1dtG = (N7/(lₙ^2))*fixed_matrices.A1*(rT*ddt)-SutdtG*GT
    
    # eqD17 in [2], matrix H2dt
    H2dtG = (N8/(lₙ^2))*fixed_matrices.A2*(rT*ddt)
    
    # -----------------------------
    # expression of C1, C2, C3, C4
    
    # eqD33 in [2], h1G
    h1G = Vec3(H1G * Ddt_e)
    
    # eqD34 in [2], h2G
    h2G = Vec3(H2G * Ddt_e)
    
    # skew-symmetric matrix of h1
    Sh1G = skew_skymmetric_matrix_from_vector(h1G)
    
    # skew-symmetric matrix of h2
    Sh2G = skew_skymmetric_matrix_from_vector(h2G)
    
    # eqD13 in [2], C1
    C1G = SWdt_e*H1G + H1dtG - H1G*Et
    
    # eqD22 in [2], C2
    C2G = SWdt_e*H2G + H2dtG - H2G*Et
    
    # eqD30 in [2], C3
    C3G = -Sh1G*GT + N7/(lₙ^2)*fixed_matrices.A1*Ddt_e*(fixed_matrices.re)' + SWdt_e*P1G*P + H1G*F1*GT

    # eqD31 in [2], C4
    C4G = -Sh2G*GT + N8/(lₙ^2)*fixed_matrices.A2*Ddt_e*(fixed_matrices.re)' + H2G*F1*GT
    
    # --------------------------------
    # expression of Udtdt, Wdt, Wdtdt
    
    # eqD12 in [2], udtdt --> Rₑ
    Udtdt = H1G*E'*ddtdt + C1G*E'*ddt
    
    # eqD20 in [2], wdtdt --> Rₑ'
    Wdt = H2G*E'*ddt
    SWdt = skew_skymmetric_matrix_from_vector(Wdt)
    
    # eqD21 in [2], wdtdt --> Rₑ'
    Wdtdt = H2G*E'*ddtdt + C2G*E'*ddt
    
    # ---------------------------------
    # expression of Tk_e, M_e and Ck_e
    ωG
    
    aux1 = (l₀/2)*ωG*conf.mat.Arho*H1G'
    wGH2GT = (l₀/2)*ωG*H2G'
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
    SIrhoeGWdt = skew_skymmetric_matrix_from_vector(IrhoeGWdt)
    term3_Ck = (aux3*IrhoeG)*H2G
    term4_Ck = -(wGH2GT*SIrhoeGWdt)*H2G
    term_viscosity_1_Ck = comp.damping*aux1*H1G
    term_viscosity_2_Ck = comp.damping*aux2*H2G
    
    # eq103 in [2]
    Ck += term1_Ck + term2_Ck + term3_Ck + term4_Ck + term_viscosity_1_Ck + term_viscosity_2_Ck
    
    # kinetic energy
    IrhoG = Rₑ*IrhoeG*Rₑ'
    wdt_G = Rₑ*Wdt
    udt_eG =  H1G*E'*ddt
    udt_G = Rₑ*udt_eG
    
    K_energy_e = K_energy_e + (l₀/2)*ωG/2*(conf.mat.Arho*(udt_G'*udt_G)) + (l₀/2)*ωG/2*(wdt_G'*(IrhoG*wdt_G))
    
    return Tk, Tdamp, M, Ck, K_energy_e, (N1, N2, utG, H1G, ωG, N7, P1G)
    
end 


# Compute beam contributions
@fastmath function get_beam_contributions!(e, allnodes, conf, sdf, fixed_matrices, comp, sol_GP, T=Float64) 
    
    # @timeit_debug to "Move to local reference system" begin
        
        # retrieve the matrix Re_0 of the beam
        Rₑ⁰ = e.R₀
        
        # information from node 1 and 2
        X₁, X₂ = local_pos(e, allnodes) 
        u₁, u₂ = local_disp(e, allnodes)
        u̇₁, u̇₂ = local_vel(e, allnodes)
        ü₁, ü₂ = local_acc(e, allnodes)
        ẇ₁, ẇ₂ = local_ang_vel(e, allnodes)
        ẅ₁, ẅ₂ = local_ang_acc(e, allnodes)
        R₁, R₂ = local_rot(e, allnodes)

        
        
        # current position of node i1 (global rs)
        x₁ =  X₁ + u₁
        x₂ =  X₂ + u₂
        
        # -------------------------------------------------------------------------------------------
        # rigidly_moving rs of the deformed configuration (v1, v2, v3)
        # -------------------------------------------------------------------------------------------
        # with respect to the new rs, the first node is placed at the origin (0) and the second at a
        # distance l: i1 --> 0 and i2 --> l
        
        # length of the beam in the undeformed configuration (CHECK)
        l₀ = e.l₀
        
        #  length of the beam in the deformed configuration
        lₙ = norm(x₂ - x₁)
        
        # eq32a in [2]: local displacement with respect to the local rs of the deformed configuration
        # (ux_0=uy_0=uz_0=uyL=uz_L=0 as a consequence of how the local rs (v1, v2, v3) is built)
        ū = lₙ - l₀

        Rₑ = local_Rₑ(x₁, x₂, Rₑ⁰[:,2])


        R̅₁ = Rₑ' * R₁ * Rₑ⁰
        R̅₂ = Rₑ' * R₂ * Rₑ⁰

        Θ̅₁ = angle_from_rotation_matrix(R̅₁)
        Θ̅₂ = angle_from_rotation_matrix(R̅₂)


        ru₁, ru₂, η, ηs, Gᵀu₁, GᵀΘ₁, Gᵀu₂, GᵀΘ₂ = auxiliary_variables(Rₑ, v₁, p, p₁, p₂)


        # inverse of the skew-symmetric matrix of Θ̅₁ and Θ̅₂
        Tₛ⁻¹Θ̅₁ = inverse_skew_skymmetric_matrix_from_angle(Θ̅₁)
        Tₛ⁻¹Θ̅₂ = inverse_skew_skymmetric_matrix_from_angle(Θ̅₂)

        P₁₁ = -Gᵀu₁ 
        P₂₁ = P₁₁
        P₁₂ = ID3-GᵀΘ₁
        P₂₂ = -GᵀΘ₁
        P₁₃ = -Gᵀu₂
        P₂₃ = P₁₃
        P₁₄ = -GᵀΘ₂
        P₂₄ = ID3-GᵀΘ₂


        # B̄′ = [r; PEᵀ]
        B̄′u₁ = ru₁
        B̄′u₂ = ru₂
        B̄′₁₁ = Rₑ * P₁₁
        B̄′₂₁ = B̄′₁₁
        B̄′₁₂ = Rₑ * P₁₂
        B̄′₂₂ = Rₑ * P₂₂
        B̄′₁₃ = Rₑ * P₁₃
        B̄′₂₃ = B̄′₁₃
        B̄′₁₄ = Rₑ * P₁₄
        B̄′₂₄ = Rₑ * P₂₄

        
        # B = B̄B̄′
        Bu₁ = B̄′u₁
        Bu₂ = B̄′u₂
        B₁₁ = Tₛ⁻¹Θ̅₁ * B̄′₁₁
        B₂₁ = Tₛ⁻¹Θ̅₂ * B̄′₂₁
        B₁₂ = Tₛ⁻¹Θ̅₁ * B̄′₁₂
        B₂₂ = Tₛ⁻¹Θ̅₂ * B̄′₂₂
        B₁₃ = Tₛ⁻¹Θ̅₁ * B̄′₁₃
        B₂₃ = Tₛ⁻¹Θ̅₂ * B̄′₂₃
        B₁₄ = Tₛ⁻¹Θ̅₁ * B̄′₁₄
        B₂₄ = Tₛ⁻¹Θ̅₂ * B̄′₂₄
        

        K̄ᵢₙₜū, K̄ᵢₙₜΘ̅ᵢᵢ, K̄ᵢₙₜΘ̅ᵢⱼ = K̄ᵢₙₜ_beam(mat, geom, l₀)

        # T̄ᵢₙₜ = K̄ᵢₙₜ D̄
        T̄ᵢₙₜū  = K̄ᵢₙₜū   * ū
        T̄ᵢₙₜΘ̅₁ = K̄ᵢₙₜΘ̅ᵢᵢ * Θ̅₁ + K̄ᵢₙₜΘ̅ᵢⱼ * Θ̅₂
        T̄ᵢₙₜΘ̅₂ = K̄ᵢₙₜΘ̅ᵢⱼ * Θ̅₁ + K̄ᵢₙₜΘ̅ᵢᵢ * Θ̅₂


        # Tᵢₙₜ = Bᵀ T̄ᵢₙₜ
        Tᵢₙₜu₁ = Bu₁*T̄ᵢₙₜū + B₁₁*T̄ᵢₙₜΘ̅₁ + B₂₁*T̄ᵢₙₜΘ̅₂
        TᵢₙₜΘ₁ =             B₁₂*T̄ᵢₙₜΘ̅₁ + B₂₂*T̄ᵢₙₜΘ̅₂
        Tᵢₙₜu₂ = Bu₂*T̄ᵢₙₜū + B₁₃*T̄ᵢₙₜΘ̅₁ + B₂₃*T̄ᵢₙₜΘ̅₂
        TᵢₙₜΘ₂ =             B₁₄*T̄ᵢₙₜΘ̅₁ + B₂₄*T̄ᵢₙₜΘ̅₂

        # [N̄ M̄′₁ M̄′₂] = B̄ T̄ᵢₙₜ
        N̄   = T̄ᵢₙₜū
        M̄′₁ = Tₛ⁻¹Θ̅₁ * T̄ᵢₙₜΘ̅₁
        M̄′₂ = Tₛ⁻¹Θ̅₂ * T̄ᵢₙₜΘ̅₂

        # Qₛ = Pᵀ [M̄′₁ M̄′₂]
        Qₛu₁ = P₁₁ * M̄′₁ + P₂₁ * M̄′₂
        QₛΘ₁ = P₁₂ * M̄′₁ + P₂₂ * M̄′₂
        Qₛu₂ = P₁₃ * M̄′₁ + P₂₃ * M̄′₂
        QₛΘ₂ = P₁₄ * M̄′₁ + P₂₄ * M̄′₂

        # Q = S(Qₛ)
        Qu₁ = skew_skymmetric_matrix_from_vector(Qₛu₁)
        QΘ₁ = skew_skymmetric_matrix_from_vector(QₛΘ₁)
        Qu₂ = skew_skymmetric_matrix_from_vector(Qₛu₂)
        QΘ₂ = skew_skymmetric_matrix_from_vector(QₛΘ₂)

        a = Vec3(0, η*(M̄′₁[1] + M̄′₂[1] - M̄′₁[2] + M̄′₂[2])/lₙ, (M̄′₁[3] + M̄′₂[3])/lₙ)

        D₃ = (ID3 - v₁*v₁')/ln

        # DN̄ (symmetric)
        DN̄₁₁ = D₃*N̄
        DN̄₃₃ = DN̄₁₁
        DN̄₁₃ = -DN̄₁₁
        DN̄₃₁ = DN̄₁₃

        # EQGᵀEᵀ (diagonal)
        EQGᵀEᵀ₁₁ = Rₑ*Qu₁*Gᵀu₁*Rₑ'
        EQGᵀEᵀ₂₂ = Rₑ*QΘ₁*GᵀΘ₁*Rₑ'
        EQGᵀEᵀ₃₃ = Rₑ*Qu₂*Gᵀu₂*Rₑ'
        EQGᵀEᵀ₄₄ = Rₑ*QΘ₂*GᵀΘ₂*Rₑ'

        # EGa (diagonal)
        # Note: Rₑ*Ga = 0 for Θ indices because Rₑ*GᵀΘ' has only non-zero values in the first column and a = [0 ...]
        EGau₁ = Rₑ*Gᵀu₁'*a
        EGau₂ = - EGau₁      # Gᵀu₂ = - Gᵀu₁

        # EGarᵀ (symmetric)
        EGarᵀ₁₁ = EGau₁*ru₁
        EGarᵀ₁₃ = -EGarᵀ₁₁  # EGau₁*ru₂ but ru₂ = - ru₁
        EGarᵀ₃₁ = EGarᵀ₁₃   # EGaru₂ = - EGaru₁ and ru₂ = - ru₁
        EGarᵀ₃₃ = -EGarᵀ₃₁

        # Kₘ = DN̄ - EQGᵀEᵀ + EGarᵀ (symmetric)
        Kₘ₁₁ = DN̄₁₁ - EQGᵀEᵀ₁₁ + EGarᵀ₁₁
        Kₘ₂₂ =      - EQGᵀEᵀ₂₂
        Kₘ₃₃ = DN̄₃₃ - EQGᵀEᵀ₃₃ + EGarᵀ₃₃
        Kₘ₄₄ =      - EQGᵀEᵀ₄₄
        Kₘ₁₃ = DN̄₁₃            + EGarᵀ₁₃
        Kₘ₃₁ = Kₘ₁₃                       # DN̄₃₁ = DN̄₁₃ and EGarᵀ₃₁ = EGarᵀ₁₃

        η₁, μ₁ = compute_η_μ(Θ̅₁)
        η₂, μ₂ = compute_η_μ(Θ̅₂)


        M̄₁ = T̄ᵢₙₜΘ̅₁
        M̄₂ = T̄ᵢₙₜΘ̅₂

        K̄ₕ₁ = compute_K̄ₕ(Θ̅₁, M̄₁, Tₛ⁻¹Θ̅₁)
        K̄ₕ₂ = compute_K̄ₕ(Θ̅₂, M̄₂, Tₛ⁻¹Θ̅₂)


        # K̃ (symmetric)

        K̃₁₁ = B̄′₁₁ * K̄ₕ₁ * B̄′₁₁  +  B̄′₂₁ * K̄ₕ₂ * B̄′₂₁  +  Kₘ₁₁
        K̃₁₂ = B̄′₁₁ * K̄ₕ₁ * B̄′₁₂  +  B̄′₂₁ * K̄ₕ₂ * B̄′₂₂
        K̃₁₃ = B̄′₁₁ * K̄ₕ₁ * B̄′₁₃  +  B̄′₂₁ * K̄ₕ₂ * B̄′₂₃  +  Kₘ₁₃
        K̃₁₄ = B̄′₁₁ * K̄ₕ₁ * B̄′₁₄  +  B̄′₂₁ * K̄ₕ₂ * B̄′₂₄

        K̃₂₁ = K̃₁₂
        K̃₂₂ = B̄′₁₂ * K̄ₕ₁ * B̄′₁₂  +  B̄′₂₂ * K̄ₕ₂ * B̄′₂₂  +  Kₘ₂₂
        K̃₂₃ = B̄′₁₂ * K̄ₕ₁ * B̄′₁₃  +  B̄′₂₂ * K̄ₕ₂ * B̄′₂₃
        K̃₂₄ = B̄′₁₂ * K̄ₕ₁ * B̄′₁₄  +  B̄′₂₂ * K̄ₕ₂ * B̄′₂₄

        K̃₃₁ = K̃₁₃
        K̃₃₂ = K̃₂₃
        K̃₃₃ = B̄′₁₃ * K̄ₕ₁ * B̄′₁₃  +  B̄′₂₃ * K̄ₕ₂ * B̄′₂₃  +  Kₘ₃₃
        K̃₃₄ = B̄′₁₃ * K̄ₕ₁ * B̄′₁₄  +  B̄′₂₃ * K̄ₕ₂ * B̄′₂₄

        K̃₄₁ = K̃₁₄
        K̃₄₂ = K̃₂₄
        K̃₄₃ = K̃₃₄
        K̃₄₄ = B̄′₁₄ * K̄ₕ₁ * B̄′₁₄  +  B̄′₂₄ * K̄ₕ₂ * B̄′₂₄  +  Kₘ₄₄

        
        # Kᵢₙₜ (symmetric)
        
        Kᵢₙₜ₁₁ = Bu₁ * K̄ᵢₙₜū * Bu₁      +      B₁₁ * K̄ᵢₙₜΘ̅ᵢᵢ * B₁₁  +  B₂₁ * K̄ᵢₙₜΘ̅ᵢⱼ * B₁₁     +      B₁₁ * K̄ᵢₙₜΘ̅ᵢⱼ * B₂₁  +  B₂₁ * K̄ᵢₙₜΘ̅ᵢᵢ * B₂₁     +      K̃₁₁
        Kᵢₙₜ₁₂ =                               B₁₁ * K̄ᵢₙₜΘ̅ᵢᵢ * B₁₂  +  B₂₁ * K̄ᵢₙₜΘ̅ᵢⱼ * B₁₂     +      B₁₁ * K̄ᵢₙₜΘ̅ᵢⱼ * B₂₂  +  B₂₁ * K̄ᵢₙₜΘ̅ᵢᵢ * B₂₂     +      K̃₁₂
        Kᵢₙₜ₁₃ = Bu₁ * K̄ᵢₙₜū * Bu₂      +      B₁₁ * K̄ᵢₙₜΘ̅ᵢᵢ * B₁₃  +  B₂₁ * K̄ᵢₙₜΘ̅ᵢⱼ * B₁₃     +      B₁₁ * K̄ᵢₙₜΘ̅ᵢⱼ * B₂₃  +  B₂₁ * K̄ᵢₙₜΘ̅ᵢᵢ * B₂₃     +      K̃₁₃
        Kᵢₙₜ₁₄ =                               B₁₁ * K̄ᵢₙₜΘ̅ᵢᵢ * B₁₄  +  B₂₁ * K̄ᵢₙₜΘ̅ᵢⱼ * B₁₄     +      B₁₁ * K̄ᵢₙₜΘ̅ᵢⱼ * B₂₄  +  B₂₁ * K̄ᵢₙₜΘ̅ᵢᵢ * B₂₄     +      K̃₁₄

        Kᵢₙₜ₂₁ = Kᵢₙₜ₁₂
        Kᵢₙₜ₂₂ =                               B₁₂ * K̄ᵢₙₜΘ̅ᵢᵢ * B₁₂  +  B₂₂ * K̄ᵢₙₜΘ̅ᵢⱼ * B₁₂     +      B₁₂ * K̄ᵢₙₜΘ̅ᵢⱼ * B₂₂  +  B₂₂ * K̄ᵢₙₜΘ̅ᵢᵢ * B₂₂     +      K̃₂₂
        Kᵢₙₜ₂₃ =                               B₁₂ * K̄ᵢₙₜΘ̅ᵢᵢ * B₁₃  +  B₂₂ * K̄ᵢₙₜΘ̅ᵢⱼ * B₁₃     +      B₁₂ * K̄ᵢₙₜΘ̅ᵢⱼ * B₂₃  +  B₂₂ * K̄ᵢₙₜΘ̅ᵢᵢ * B₂₃     +      K̃₂₃
        Kᵢₙₜ₂₄ =                               B₁₂ * K̄ᵢₙₜΘ̅ᵢᵢ * B₁₄  +  B₂₂ * K̄ᵢₙₜΘ̅ᵢⱼ * B₁₄     +      B₁₂ * K̄ᵢₙₜΘ̅ᵢⱼ * B₂₄  +  B₂₂ * K̄ᵢₙₜΘ̅ᵢᵢ * B₂₄     +      K̃₂₄

        Kᵢₙₜ₃₁ = Kᵢₙₜ₁₃
        Kᵢₙₜ₃₂ = Kᵢₙₜ₂₃
        Kᵢₙₜ₃₃ = Bu₂ * K̄ᵢₙₜū * Bu₂      +      B₁₃ * K̄ᵢₙₜΘ̅ᵢᵢ * B₁₃  +  B₂₃ * K̄ᵢₙₜΘ̅ᵢⱼ * B₁₃     +      B₁₃ * K̄ᵢₙₜΘ̅ᵢⱼ * B₂₃  +  B₂₃ * K̄ᵢₙₜΘ̅ᵢᵢ * B₂₃     +      K̃₃₃
        Kᵢₙₜ₃₄ =                               B₁₃ * K̄ᵢₙₜΘ̅ᵢᵢ * B₁₄  +  B₂₃ * K̄ᵢₙₜΘ̅ᵢⱼ * B₁₄     +      B₁₃ * K̄ᵢₙₜΘ̅ᵢⱼ * B₂₄  +  B₂₃ * K̄ᵢₙₜΘ̅ᵢᵢ * B₂₄     +      K̃₃₄

        Kᵢₙₜ₄₁ = Kᵢₙₜ₁₄
        Kᵢₙₜ₄₂ = Kᵢₙₜ₂₄
        Kᵢₙₜ₄₃ = Kᵢₙₜ₃₄
        Kᵢₙₜ₄₄ =                               B₁₄ * K̄ᵢₙₜΘ̅ᵢᵢ * B₁₄  +  B₂₄ * K̄ᵢₙₜΘ̅ᵢⱼ * B₁₄     +      B₁₄ * K̄ᵢₙₜΘ̅ᵢⱼ * B₂₄  +  B₂₄ * K̄ᵢₙₜΘ̅ᵢᵢ * B₂₄     +      K̃₄₄
        

        
        # Internal energy    
        Phi_energy_e = (ū*T̄ᵢₙₜū + dot(Θ̅₁, T̄ᵢₙₜΘ̅₁) + dot(Θ̅₂, T̄ᵢₙₜΘ̅₂))/2
        

        

        

        U̇₁ = Rₑ' * u̇₁
        U̇₂ = Rₑ' * u̇₂
        Ẇ₁ = Rₑ' * ẇ₁
        Ẇ₂ = Rₑ' * ẇ₂

        SU̇₁ = skew_skymmetric_matrix_from_vector(U̇₁)
        SU̇₂ = skew_skymmetric_matrix_from_vector(U̇₂)
        SẆ₁ = skew_skymmetric_matrix_from_vector(Ẇ₁)
        SẆ₂ = skew_skymmetric_matrix_from_vector(Ẇ₂)

        ẆᵉU̇₁ = Gᵀu₁*U̇₁
        ẆᵉU̇₂ = Gᵀu₂*U̇₂
        ẆᵉẆ₁ = GᵀΘ₁*Ẇ₁
        ẆᵉẆ₂ = GᵀΘ₁*Ẇ₂
        
        SẆᵉU̇₁ = skew_skymmetric_matrix_from_vector(ẆᵉU̇₁)
        SẆᵉU̇₂ = skew_skymmetric_matrix_from_vector(ẆᵉU̇₂)
        SẆᵉẆ₁ = skew_skymmetric_matrix_from_vector(ẆᵉẆ₁)
        SẆᵉẆ₂ = skew_skymmetric_matrix_from_vector(ẆᵉẆ₂)
        
        
        # initialise the local matrices used in the Gauss loop
        K_energy_e = zero(T)
        C_energy_e = zero(T)
        
        Tₖ = zeros(Vec12{T})
        Tᵥ = zeros(Vec12{T})
        M = zeros(Mat1212{T}) 
        Cₖ = zeros(Mat1212{T})
        Tc = zeros(Vec12{T}) 
        Kc = zeros(Mat1212{T}) 
        Cc = zeros(Mat1212{T}) 

        # cycle among the Gauss positions
        for iG in 1:comp.nG

            zG = comp.zG[iG]
            ωG = comp.ωG[iG]

            ξ = l₀*(zG+1)/2
    
            # Shape functions
            N₁ = 1-ξ/l₀
            N₂ = 1-N₁
            N₃ = ξ*(1-ξ/l₀)^2
            N₄ = -(1-ξ/l₀)*((ξ^2)/l₀)
            N₅ = (1-3*ξ/l₀)*(1-ξ/l₀)
            N₆ = (3*ξ/l₀-2)*(ξ/l₀)
            N₇ = N3+N4
            N₈ = N5+N6-1


            uᵗ = Vec3(0, N₃*Θ̄₁[3] + N₄*Θ̄₂[3], -N₃*Θ̄₁[2] + -N₄*Θ̄₂[2])
            Θ̄  = Vec3(N₁*Θ̄₁[1] + N₂*Θ̄₂[1], N₅*Θ̄₁[2] + N₆*Θ̄₂[2], N₅*Θ̄₁[3] + N₆*Θ̄₂[3])

            Suᵗ = skew_skymmetric_matrix_from_vector(uᵗ)
            SΘ̄ = skew_skymmetric_matrix_from_vector(Θ̄)

            R̄ = ID3 + SΘ̄

            Īᵨ = R̄*conf.mat.Jᵨ*R̄'

            P₁Pu₁, P₁PΘ₁, P₁Pu₂, P₁PΘ₂, P₂Pu₁, P₂PΘ₁, P₂Pu₂, P₂PΘ₂ = Pmatrices(N₁, N₂, N₃, N₄, N₅, N₆, lₙ, η, ηs...)

            H₁u₁ = N₁*ID3 + P₁Pu₁ - Suᵗ*Gᵀu₁
            H₁Θ₁ =          P₁PΘ₁ - Suᵗ*GᵀΘ₁
            H₁u₂ = N₂*ID3 + P₁Pu₂ - Suᵗ*Gᵀu₂
            H₁Θ₂ =          P₁PΘ₂ - Suᵗ*GᵀΘ₂

            H₂u₁ = P₂Pu₁ + Gᵀu₁
            H₂Θ₁ = P₂PΘ₁ + GᵀΘ₁
            H₂u₂ = P₂Pu₂ + Gᵀu₂
            H₂Θ₂ = P₂PΘ₂ + GᵀΘ₂

            u̇ᵗu₁ =  P₁Pu₁ * U̇₁
            u̇ᵗΘ₁ =  P₁Pu₁ * Ẇ₁
            u̇ᵗu₂ =  P₁Pu₂ * U̇₂
            u̇ᵗΘ₂ =  P₁Pu₂ * Ẇ₂

            rḋ = dot(ru₁, u̇₁) + dot(ru₂, u̇₂)
            
            N₇rḋ = N₇/lₙ^2 * rḋ
            Ḣ₁u₁ = N₇rḋ * Diagonal(Vec3(0, -1, -1)) - skew_skymmetric_matrix_from_vector(u̇ᵗu₁) * Gᵀu₁
            Ḣ₁Θ₁ =                                  - skew_skymmetric_matrix_from_vector(u̇ᵗΘ₁) * GᵀΘ₁
            Ḣ₁u₂ = N₇rḋ * Diagonal(Vec3(0,  1,  1)) - skew_skymmetric_matrix_from_vector(u̇ᵗu₂) * Gᵀu₂
            Ḣ₁Θ₂ =                                  - skew_skymmetric_matrix_from_vector(u̇ᵗΘ₂) * GᵀΘ₂

            N₈rḋ = N₈/lₙ^2 * rḋ
            Ḣ₂u₁ = @SMatrix [0 0 0; 0 0 N₈rḋ; 0 -N₈rḋ 0]
            Ḣ₂u₁ = -Ḣ₂u₁

            h₁ = H₁u₁ * U̇₁ + H₁Θ₁ * Ẇ₁ + H₁u₂ * U̇₂ + H₁Θ₂ * Ẇ₂
            h₂ = H₂u₁ * U̇₁ + H₂Θ₁ * Ẇ₁ + H₂u₂ * U̇₂ + H₂Θ₂ * Ẇ₂
            Sh₁ = skew_skymmetric_matrix_from_vector(h₁)
            Sh₂ = skew_skymmetric_matrix_from_vector(h₂)

            Ẇᵉ = Gᵀu₁ * U̇₁ + GᵀΘ₁ * Ẇ₁ + Gᵀu₂ * U̇₂ + GᵀΘ₂ * Ẇ₂
            SẆᵉ = skew_skymmetric_matrix_from_vector(Ẇᵉ)

            C₁u₁ = SẆᵉ * H₁u₁ + Ḣ₁u₁ - H₁u₁ * SẆᵉ
            C₁Θ₁ = SẆᵉ * H₁Θ₁ + Ḣ₁Θ₁ - H₁Θ₁ * SẆᵉ
            C₁u₂ = SẆᵉ * H₁u₂ + Ḣ₁u₂ - H₁u₂ * SẆᵉ
            C₁Θ₂ = SẆᵉ * H₁Θ₂ + Ḣ₁Θ₂ - H₁Θ₂ * SẆᵉ

            C₂u₁ = SẆᵉ * H₂u₁ + Ḣ₁u₁ - H₂u₁ * SẆᵉ
            C₂Θ₁ = SẆᵉ * H₂Θ₁ + Ḣ₁Θ₁ - H₂Θ₁ * SẆᵉ
            C₂u₂ = SẆᵉ * H₂u₂ + Ḣ₁u₂ - H₂u₂ * SẆᵉ
            C₂Θ₂ = SẆᵉ * H₂Θ₂ + Ḣ₁Θ₂ - H₂Θ₂ * SẆᵉ

            H₁F₁ = H₁u₁ * SU̇₁ + H₁Θ₁ * SẆ₁ + H₁u₂ * SU̇₂ + H₁Θ₂ * SẆ₂
            H₂F₁ = H₂u₁ * SU̇₁ + H₂Θ₁ * SẆ₁ + H₂u₂ * SU̇₂ + H₂Θ₂ * SẆ₂

            A₁ḊrEu₁ = @SMatrix [0 0 0; U̇₁[2]-U̇₂[2] 0 0; U̇₁[3]-U̇₂[3] 0 0]
            A₁ḊrEu₂ = @SMatrix [0 0 0; -U̇₁[2]+U̇₂[2] 0 0; -U̇₁[3]+U̇₂[3] 0 0]
            C₃u₁ = -Sh₁*Gᵀu₁ + N₇/lₙ^2*A₁ḊrEu₁ + SẆᵉ*P₁Pu₁ + H₁F₁*Gᵀu₁
            C₃Θ₁ = -Sh₁*GᵀΘ₁ +                 SẆᵉ*P₁PΘ₁ + H₁F₁*GᵀΘ₁
            C₃u₂ = -Sh₁*Gᵀu₂ + N₇/lₙ^2*A₁ḊrEu₂ + SẆᵉ*P₁Pu₂ + H₁F₁*Gᵀu₂
            C₃Θ₂ = -Sh₁*GᵀΘ₂ +                 SẆᵉ*P₁PΘ₂ + H₁F₁*GᵀΘ₂

            A₂ḊrEu₁ = @SMatrix [0 0 0; -U̇₁[3]+U̇₂[3] 0 0; U̇₁[2]-U̇₂[2] 0 0]
            A₂ḊrEu₂ = @SMatrix [0 0 0; U̇₁[3]-U̇₂[3] 0 0; -U̇₁[2]+U̇₂[2] 0 0]
            C₄u₁ = -Sh₂*Gᵀu₁ + N₈/lₙ^2*A₂ḊrEu₁ + H₂F₁*Gᵀu₁
            C₄Θ₁ = -Sh₂*GᵀΘ₁                   + H₂F₁*GᵀΘ₁
            C₄u₂ = -Sh₂*Gᵀu₂ + N₈/lₙ^2*A₂ḊrEu₂ + H₂F₁*Gᵀu₂
            C₄Θ₂ = -Sh₂*GᵀΘ₂                   + H₂F₁*GᵀΘ₂

            H₁Eᵀḋ = H₁u₁ * Rₑ' * U̇₁ + H₁Θ₁ * Rₑ' * Ẇ₁ + H₁u₂ * Rₑ' * U̇₂ + H₁Θ₂ * Rₑ' * Ẇ₂
            u̇₀ = Rₑ * H₁Eᵀḋ

            H₁Eᵀd̈ = H₁u₁ * Rₑ' * Ü₁ + H₁Θ₁ * Rₑ' * Ẅ₁ + H₁u₂ * Rₑ' * Ü₂ + H₁Θ₂ * Rₑ' * Ẅ₂
            C₁Eᵀḋ = C₁u₁ * Rₑ' * Ü₁ + C₁Θ₁ * Rₑ' * Ẅ₁ + C₁u₂ * Rₑ' * Ü₂ + C₁Θ₂ * Rₑ' * Ẅ₂
            ü₀ = Rₑ * H₁Eᵀd̈ + C₁Eᵀḋ

            H₂Eᵀḋ = H₂u₁ * Rₑ' * U̇₁ + H₂Θ₁ * Rₑ' * Ẇ₁ + H₂u₂ * Rₑ' * U̇₂ + H₂Θ₂ * Rₑ' * Ẇ₂
            ẇ₀ = Rₑ * H₂Eᵀḋ

            H₂Eᵀd̈ = H₂u₁ * Rₑ' * Ü₁ + H₂Θ₁ * Rₑ' * Ẅ₁ + H₂u₂ * Rₑ' * Ü₂ + H₂Θ₂ * Rₑ' * Ẅ₂
            C₂Eᵀḋ = C₂u₁ * Rₑ' * Ü₁ + C₂Θ₁ * Rₑ' * Ẅ₁ + C₂u₂ * Rₑ' * Ü₂ + C₂Θ₂ * Rₑ' * Ẅ₂
            ẅ₀ = Rₑ * H₂Eᵀd̈ + C₂Eᵀḋ

            Ẇ₀ = H₂Eᵀḋ
            SẆ₀ = skew_skymmetric_matrix_from_vector(Ẇ₀)
            Rₑᵀü₀ = Rₑ'*ü₀
            ĪᵨRₑᵀẅ₀ = Īᵨ*Rₑ'*ẅ₀
            SẆ₀ĪᵨẆ₀ = SẆ₀*Īᵨ*Ẇ₀
            Tₖu₁ += ω*l₀/2 * (Aᵨ*H₁u₁'*Rₑᵀü₀ + H₂u₁'*ĪᵨRₑᵀẅ₀ + H₂u₁'*SẆ₀ĪᵨẆ₀)
            TₖΘ₁ += ω*l₀/2 * (Aᵨ*H₁Θ₁'*Rₑᵀü₀ + H₂Θ₁'*ĪᵨRₑᵀẅ₀ + H₂Θ₁'*SẆ₀ĪᵨẆ₀)
            Tₖu₂ += ω*l₀/2 * (Aᵨ*H₁u₂'*Rₑᵀü₀ + H₂u₂'*ĪᵨRₑᵀẅ₀ + H₂u₂'*SẆ₀ĪᵨẆ₀)
            TₖΘ₂ += ω*l₀/2 * (Aᵨ*H₁Θ₂'*Rₑᵀü₀ + H₂Θ₂'*ĪᵨRₑᵀẅ₀ + H₂Θ₂'*SẆ₀ĪᵨẆ₀)


            M₁₁ += ω*l₀/2 * (Aᵨ*H₁u₁'*H₁u₁ + H₂u₁'*Īᵨ*H₂u₁)
            M₁₂ += ω*l₀/2 * (Aᵨ*H₁u₁'*H₁Θ₁ + H₂u₁'*Īᵨ*H₂Θ₁)
            M₁₃ += ω*l₀/2 * (Aᵨ*H₁u₁'*H₁u₂ + H₂u₁'*Īᵨ*H₂u₂)
            M₁₄ += ω*l₀/2 * (Aᵨ*H₁u₁'*H₁Θ₂ + H₂u₁'*Īᵨ*H₂Θ₂)

            M₂₁ += M₁₂
            M₂₂ += ω*l₀/2 * (Aᵨ*H₁Θ₁'*H₁Θ₁ + H₂Θ₁'*Īᵨ*H₂Θ₁)
            M₂₃ += ω*l₀/2 * (Aᵨ*H₁Θ₁'*H₁u₂ + H₂Θ₁'*Īᵨ*H₂u₂)
            M₂₄ += ω*l₀/2 * (Aᵨ*H₁Θ₁'*H₁Θ₂ + H₂Θ₁'*Īᵨ*H₂Θ₂)

            M₃₁ += M₁₃
            M₃₂ += M₂₃
            M₃₃ += ω*l₀/2 * (Aᵨ*H₁u₂'*H₁u₂ + H₂u₂'*Īᵨ*H₂u₂)
            M₃₄ += ω*l₀/2 * (Aᵨ*H₁u₂'*H₁Θ₂ + H₂u₂'*Īᵨ*H₂Θ₂)

            M₄₁ += M₁₄
            M₄₂ += M₂₄
            M₄₃ += M₃₄
            M₄₄ += ω*l₀/2 * (Aᵨ*H₁Θ₂'*H₁Θ₂ + H₂Θ₂'*Īᵨ*H₂Θ₂)









            # eqD4 in [2]: local transerve displacement vector
        utG = reshapeP1G(N3, N4, Θ̅₁, Θ̅₂)
        SutG = skew_skymmetric_matrix_from_vector(utG)
        
        # eqD14 in [2]
        utdtG = P1G*P * E' * ddt
        SutdtG = skew_skymmetric_matrix_from_vector(utdtG)
        
        # ???: local transerve rotational vector
        ΘG_bar = P2G * Vec6(Θ̅₁[1], Θ̅₁[2], Θ̅₁[3],  Θ̅₂[1], Θ̅₂[2], Θ̅₂[3])
        SuΘG_bar = skew_skymmetric_matrix_from_vector(ΘG_bar)
        
        # eq74 in [2]: RG_bar, rotation matrix associated with the local cross-section rotation
        RG_bar = ID3 + SuΘG_bar
        
        # eq100 in [2]
        IrhoeG = RG_bar*(conf.mat.Jᵨ)*RG_bar'
        
        # ------------------------------------
        # expression of H1, H2 and H1dt, H2dt
        
        # eqD6 in [2], matrix H1
        H1G = NG + P1G*P - SutG*GT
        
        # eqD16 in [2], matrix H2
        H2G = P2G*P + GT
        
        # eqD14 in [2], matrix H1dt
        H1dtG = (N7/(lₙ^2))*fixed_matrices.A1*(rT*ddt)-SutdtG*GT
        
        # eqD17 in [2], matrix H2dt
        H2dtG = (N8/(lₙ^2))*fixed_matrices.A2*(rT*ddt)
        
        # -----------------------------
        # expression of C1, C2, C3, C4
        
        # eqD33 in [2], h1G
        h1G = Vec3(H1G * Ddt_e)
        
        # eqD34 in [2], h2G
        h2G = Vec3(H2G * Ddt_e)
        
        # skew-symmetric matrix of h1
        Sh1G = skew_skymmetric_matrix_from_vector(h1G)
        
        # skew-symmetric matrix of h2
        Sh2G = skew_skymmetric_matrix_from_vector(h2G)
        
        # eqD13 in [2], C1
        C1G = SWdt_e*H1G + H1dtG - H1G*Et
        
        # eqD22 in [2], C2
        C2G = SWdt_e*H2G + H2dtG - H2G*Et
        
        # eqD30 in [2], C3
        C3G = -Sh1G*GT + N7/(lₙ^2)*fixed_matrices.A1*Ddt_e*(fixed_matrices.re)' + SWdt_e*P1G*P + H1G*F1*GT

        # eqD31 in [2], C4
        C4G = -Sh2G*GT + N8/(lₙ^2)*fixed_matrices.A2*Ddt_e*(fixed_matrices.re)' + H2G*F1*GT
        
        # --------------------------------
        # expression of Udtdt, Wdt, Wdtdt
        
        # eqD12 in [2], udtdt --> Rₑ
        Udtdt = H1G*E'*ddtdt + C1G*E'*ddt
        
        # eqD20 in [2], wdtdt --> Rₑ'
        Wdt = H2G*E'*ddt
        SWdt = skew_skymmetric_matrix_from_vector(Wdt)
        
        # eqD21 in [2], wdtdt --> Rₑ'
        Wdtdt = H2G*E'*ddtdt + C2G*E'*ddt
        
        # ---------------------------------
        # expression of Tk_e, M_e and Ck_e
        ωG
        
        aux1 = (l₀/2)*ωG*conf.mat.Arho*H1G'
        wGH2GT = (l₀/2)*ωG*H2G'
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
        SIrhoeGWdt = skew_skymmetric_matrix_from_vector(IrhoeGWdt)
        term3_Ck = (aux3*IrhoeG)*H2G
        term4_Ck = -(wGH2GT*SIrhoeGWdt)*H2G
        term_viscosity_1_Ck = comp.damping*aux1*H1G
        term_viscosity_2_Ck = comp.damping*aux2*H2G
        
        # eq103 in [2]
        Ck += term1_Ck + term2_Ck + term3_Ck + term4_Ck + term_viscosity_1_Ck + term_viscosity_2_Ck
        
        # kinetic energy
        IrhoG = Rₑ*IrhoeG*Rₑ'
        wdt_G = Rₑ*Wdt
        udt_eG =  H1G*E'*ddt
        udt_G = Rₑ*udt_eG
        
        K_energy_e = K_energy_e + (l₀/2)*ωG/2*(conf.mat.Arho*(udt_G'*udt_G)) + (l₀/2)*ωG/2*(wdt_G'*(IrhoG*wdt_G))
        
        return Tk, Tdamp, M, Ck, K_energy_e, (N1, N2, utG, H1G, ωG, N7, P1G)

            # @timeit_debug to "Dynamic contributions" begin

                Tk, Tdamp, M, Ck, K_energy_e, contact_vals = compute_dynamic_contribution(zG_v[iG], l₀,  wG_v[iG], fixed_matrices.P1G_v[iG], fixed_matrices.P2G_v[iG], fixed_matrices.NG_v[iG], Θ̅₁, Θ̅₂, P, E, ddt, GT, lₙ, rT, Ddt_e, SWdt_e, Et, F1, ddtdt, Rₑ, Tk, Tdamp, M, Ck, K_energy_e, fixed_matrices, conf, comp)
            # end 

            # @timeit_debug to "Contact contributions" begin


                if !isnothing(sdf)
                    Tct, Kct, Cc, C_energy_e =  compute_contact_contribution(contact_vals, x₁, x₂, Rₑ, comp, sdf, E, ddt, l₀, EG, lₙ, fixed_matrices, rT, P, Tct, Kct, Cc, C_energy_e, iG, e, sol_GP, T)
                end     

            # end 
            
        end 

    # end 
         
    # @timeit_debug to "Return to global reference system" begin
        
        Tk = E * Tk
        Tdamp =  E * Tdamp
        M = E * M * E'
        Ck = E * Ck * E'
        

        # Note: this has no effect on convergence rate from my tests
        
        # Deltn1_1, Deltn1_2 = local_rot_delt(e, allnodes)
        
        # thi_n1_1 = angle_from_rotation_matrix(Deltn1_1)
        # thi_n1_2 = angle_from_rotation_matrix(Deltn1_2)
        
        # Tₛ⁻¹_th1g = inverse_skew_skymmetric_matrix_from_angle(thi_n1_1)
        # Tₛ⁻¹_th2g = inverse_skew_skymmetric_matrix_from_angle(thi_n1_2)
        
        # Bt = compute_Bt(Tₛ⁻¹_th1g, Tₛ⁻¹_th2g)
        
        # M = M*Bt
        # Ck = Ck*Bt
        # Cc = Cc*Bt
        
        Ck += Cc 
        
    # end 

    
    return Kint, Tint, Tk, Tdamp, Ck, M, Tct, Kct, Phi_energy_e, K_energy_e, C_energy_e
    
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

        @timeit_debug to "Compute elemental contributions" begin

        Kint_buf = [zeros(Mat1212) for _ in 1:length(allbeams)]
        Ck_buf = [zeros(Mat1212) for _ in 1:length(allbeams)]
        M_buf = [zeros(Mat1212) for _ in 1:length(allbeams)]
        Kct_buf = [zeros(Mat1212) for _ in 1:length(allbeams)]
        Tint_buf = [zeros(Vec12) for _ in 1:length(allbeams)]
        Tk_buf = [zeros(Vec12) for _ in 1:length(allbeams)]
        Tdamp_buf = [zeros(Vec12) for _ in 1:length(allbeams)]
        Tct_buf = [zeros(Vec12) for _ in 1:length(allbeams)]
        Phi_energy_buf = [zero(T) for _ in 1:length(allbeams)]
        K_energy_buf = [zero(T) for _ in 1:length(allbeams)]
        C_energy_buf = [zero(T) for _ in 1:length(allbeams)]

        for i in 1:length(allbeams)
        
            #----------------------------------------
            # Compute the contibution from the e beam
            e = allbeams[i]
            #@infiltrate
            Kint_buf[i], Tint_buf[i], Tk_buf[i], Tdamp_buf[i], Ck_buf[i], M_buf[i], Tct_buf[i], Kct_buf[i], Phi_energy_buf[i], K_energy_buf[e.ind], C_energy_buf[e.ind]  = get_beam_contributions!(e, allnodes, conf, sdf, fixed_matrices, comp, sol_GP, T) 

        end

        end

        @timeit_debug to "Assemble elemental contributions" begin

        for (i, e) in enumerate(LazyRows(allbeams))

            #-----------------------
            # Assemble contributions
            
            i1 = e.node1
            i2 = e.node2
            idof1 = allnodes.idof_6[i1]
            idof2 = allnodes.idof_6[i2]
            
            idof = vcat(idof1, idof2)
            
            energy.Phi_energy +=  Phi_energy_buf[i]
            energy.K_energy += K_energy_buf[i]
            energy.C_energy +=  C_energy_buf[i]
        
            update_spmat_sum(matrices.Kint, e.sparsity_map, Kint_buf[i])
            update_spmat_sum(matrices.Ck, e.sparsity_map, Ck_buf[i])
            update_spmat_sum(matrices.M, e.sparsity_map, M_buf[i])
            update_spmat_sum(matrices.Kct, e.sparsity_map, Kct_buf[i])

            
            update_vec_sum(matrices.Tk, idof, Tk_buf[i])
            update_vec_sum(matrices.Tdamp, idof, Tdamp_buf[i])
            update_vec_sum(matrices.Tint, idof, Tint_buf[i])
            update_vec_sum(matrices.Tct, idof, Tct_buf[i])

        end
                               
        end



    end 


    
end 