function get_beam_info(beam, nodes)
     
    Nₐ = beam.node1
    Nᵦ = beam.node2

    # information from node 1 and 2
    X₁, X₂ = nodes.X₀[Nₐ], nodes.X₀[Nᵦ]
    u₁, u₂ = nodes.u[Nₐ], nodes.u[Nᵦ]
    R₁, R₂ = nodes.R[Nₐ], nodes.R[Nᵦ]

    init = (X₁, X₂, beam.l₀, beam.Rₑ⁰)      
    beaminfo = (u₁, u₂, R₁, R₂, init)

    return beaminfo

end 

function beam2beam_get_RₑH₁Rₑᵀ_and_RₑdH₁Rₑᵀ(ξ, u₁, u₂, R₁, R₂, init) where T

    X₁, X₂, l₀, Rₑ⁰ = init
    
    x₁ =  X₁ + u₁
    x₂ =  X₂ + u₂
    
    lₙ = norm(x₂-x₁)

    Rₑ, _, _, _, Gᵀ¹, Gᵀ², Gᵀ³, Gᵀ⁴, _ = local_Rₑ_and_aux(x₁, x₂, R₁, R₂, Rₑ⁰[:,2], lₙ)

    R̅₁ = Rₑ' * R₁ * Rₑ⁰
    R̅₂ = Rₑ' * R₂ * Rₑ⁰

    Θ̅₁ = toangle(R̅₁)
    Θ̅₂ = toangle(R̅₂)

    N₁ = 1-ξ/l₀
    N₂ = 1-N₁
    N₃ = ξ*(1-ξ/l₀)^2
    N₄ = -(1-ξ/l₀)*((ξ^2)/l₀)

    dN₁ = -1/l₀
    dN₂ = 1/l₀
    dN₃ = (1-ξ/l₀)^2 + (-1/l₀)*2*ξ*(1-ξ/l₀)
    dN₄ = (1/l₀)*((ξ^2)/l₀) - (1-ξ/l₀)*2*ξ/l₀

    uᵗ = @SVector [0, N₃*Θ̅₁[3] + N₄*Θ̅₂[3], -N₃*Θ̅₁[2] + -N₄*Θ̅₂[2]]
    duᵗ = @SVector [0, dN₃*Θ̅₁[3] + dN₄*Θ̅₂[3], -dN₃*Θ̅₁[2] + -dN₄*Θ̅₂[2]]

    Suᵗ = skew(uᵗ)
    Sduᵗ = skew(duᵗ)

    P₁P¹ = @SMatrix [0 0 0; 0 (N₃+N₄)/lₙ 0;0 0 (N₃+N₄)/lₙ]
    P₁P² = @SMatrix [0 0 0; 0 0 N₃;0 -N₃ 0]
    P₁P³ = -P₁P¹
    P₁P⁴ = @SMatrix [0 0 0; 0 0 N₄;0 -N₄ 0]

    H₁¹ = N₁*ID3 + P₁P¹ - Suᵗ*Gᵀ¹
    H₁² =          P₁P² - Suᵗ*Gᵀ²
    H₁³ = N₂*ID3 + P₁P³ - Suᵗ*Gᵀ³
    H₁⁴ =          P₁P⁴ - Suᵗ*Gᵀ⁴

    dP₁P¹ = @SMatrix [0 0 0; 0 (dN₃+dN₄)/lₙ 0;0 0 (dN₃+dN₄)/lₙ]
    dP₁P² = @SMatrix [0 0 0; 0 0 dN₃;0 -dN₃ 0]
    dP₁P³ = -P₁P¹
    dP₁P⁴ = @SMatrix [0 0 0; 0 0 dN₄;0 -dN₄ 0]

    dH₁¹ = dN₁*ID3 + dP₁P¹ - Sduᵗ*Gᵀ¹
    dH₁² =           dP₁P² - Sduᵗ*Gᵀ²
    dH₁³ = dN₂*ID3 + dP₁P³ - Sduᵗ*Gᵀ³
    dH₁⁴ =           dP₁P⁴ - Sduᵗ*Gᵀ⁴

    RₑH₁¹Rₑᵀ = Rₑ * H₁¹ * Rₑ'
    RₑH₁²Rₑᵀ = Rₑ * H₁² * Rₑ'
    RₑH₁³Rₑᵀ = Rₑ * H₁³ * Rₑ'
    RₑH₁⁴Rₑᵀ = Rₑ * H₁⁴ * Rₑ'

    RₑdH₁¹Rₑᵀ = Rₑ * dH₁¹ * Rₑ'
    RₑdH₁²Rₑᵀ = Rₑ * dH₁² * Rₑ'
    RₑdH₁³Rₑᵀ = Rₑ * dH₁³ * Rₑ'
    RₑdH₁⁴Rₑᵀ = Rₑ * dH₁⁴ * Rₑ'

    return RₑH₁¹Rₑᵀ, RₑH₁²Rₑᵀ, RₑH₁³Rₑᵀ, RₑH₁⁴Rₑᵀ, RₑdH₁¹Rₑᵀ, RₑdH₁²Rₑᵀ, RₑdH₁³Rₑᵀ, RₑdH₁⁴Rₑᵀ

end 

function beam2beam_get_derivatives(ξ, u₁, u₂, R₁, R₂, init) where T

    X₁, X₂, l₀, Rₑ⁰ = init
    
    x₁ =  X₁ + u₁
    x₂ =  X₂ + u₂
    
    lₙ = norm(x₂ - x₁)

    Rₑ, _, _, _, _, _, _, _, _ = local_Rₑ_and_aux(x₁, x₂, R₁, R₂, Rₑ⁰[:,2], lₙ)

    R̅₁ = Rₑ' * R₁ * Rₑ⁰
    R̅₂ = Rₑ' * R₂ * Rₑ⁰

    Θ̅₁ = toangle(R̅₁)
    Θ̅₂ = toangle(R̅₂)

    N₁ = 1-ξ/l₀
    N₂ = 1-N₁
    N₃ = ξ*(1-ξ/l₀)^2
    N₄ = -(1-ξ/l₀)*((ξ^2)/l₀)

    dN₃ = (1-ξ/l₀)^2 + (-1/l₀)*2*ξ*(1-ξ/l₀)
    dN₄ = (1/l₀)*((ξ^2)/l₀) - (1-ξ/l₀)*2*ξ/l₀

    ddN₃ = 2*(1-ξ/l₀)*(-1/l₀) - 1/l₀*2*(1-ξ/l₀) - 1/l₀*2*ξ*(-1/l₀)
    ddN₄ = (1/l₀)*(1/l₀)*2*ξ + 1/l₀*2*ξ/l₀ - (1-ξ/l₀)*2/l₀

    uᵗ = @SVector [0, N₃*Θ̅₁[3] + N₄*Θ̅₂[3], -N₃*Θ̅₁[2] + -N₄*Θ̅₂[2]]
    duᵗ = @SVector [0, dN₃*Θ̅₁[3] + dN₄*Θ̅₂[3], -dN₃*Θ̅₁[2] + -dN₄*Θ̅₂[2]]
    dduᵗ = @SVector [0, ddN₃*Θ̅₁[3] + ddN₄*Θ̅₂[3], -ddN₃*Θ̅₁[2] + -ddN₄*Θ̅₂[2]]

    xᴳ = N₁*x₁ + N₂*x₂ + Rₑ*uᵗ
    dxᴳ =  (1/l₀)*(x₂-x₁) + Rₑ*duᵗ
    ddxᴳ = Rₑ*dduᵗ

    return xᴳ, dxᴳ, ddxᴳ

end 

function beam2beam_compute_dist_candidate(beaminfoₐ, beaminfoᵦ)

    ξᶜₐ = 0
    ξᶜᵦ = 0

    (u₁ₐ, u₂ₐ, R₁ₐ, R₂ₐ, initₐ) = beaminfoₐ
    (u₁ᵦ, u₂ᵦ, R₁ᵦ, R₂ᵦ, initᵦ) = beaminfoᵦ

    _, _, l₀ₐ, _ = initₐ
    _, _, l₀ᵦ, _ = initᵦ
    
    xᶜₐ, dxᶜₐ, ddxᶜₐ  = beam2beam_get_derivatives(ξᶜₐ, u₁ₐ, u₂ₐ, R₁ₐ, R₂ₐ, initₐ)
    xᶜᵦ, dxᶜᵦ, ddxᶜᵦ  = beam2beam_get_derivatives(ξᶜᵦ, u₁ᵦ, u₂ᵦ, R₁ᵦ, R₂ᵦ, initᵦ)

    A = Mat22(dot(dxᶜₐ, dxᶜₐ), dot(dxᶜₐ,dxᶜᵦ), -dot(dxᶜₐ,dxᶜᵦ), -dot(dxᶜᵦ,dxᶜᵦ))
    b = Vec2(-dot((xᶜₐ-xᶜᵦ),dxᶜₐ), -dot((xᶜₐ-xᶜᵦ),dxᶜᵦ))

    ξᶜ = A\b

    ξᶜₐ = ξᶜ[1]
    ξᶜᵦ = ξᶜ[2]

    tol = 1E-6
    norm_res = 1E6
    max_it = 50
    k = 0

    while norm_res>tol && k<max_it

        xᶜₐ, dxᶜₐ, ddxᶜₐ  = beam2beam_get_derivatives(ξᶜₐ, u₁ₐ, u₂ₐ, R₁ₐ, R₂ₐ, initₐ)
        xᶜᵦ, dxᶜᵦ, ddxᶜᵦ  = beam2beam_get_derivatives(ξᶜᵦ, u₁ᵦ, u₂ᵦ, R₁ᵦ, R₂ᵦ, initᵦ)
        
        A = Mat22(dot(dxᶜₐ, dxᶜₐ)+dot((xᶜₐ-xᶜᵦ),ddxᶜₐ), dot(dxᶜₐ,dxᶜᵦ), -dot(dxᶜₐ,dxᶜᵦ), -dot(dxᶜᵦ,dxᶜᵦ) + dot((xᶜₐ-xᶜᵦ),ddxᶜᵦ))
        b = Vec2(-dot((xᶜₐ-xᶜᵦ),dxᶜₐ), -dot((xᶜₐ-xᶜᵦ),dxᶜᵦ))
                
        Δξᶜ = A\b

        Δξᶜₐ = Δξᶜ[1]
        Δξᶜᵦ = Δξᶜ[2]

        ξᶜₐ = ξᶜₐ + Δξᶜₐ
        ξᶜᵦ = ξᶜᵦ + Δξᶜᵦ
        
        norm_res = norm(b)
        k = k+1

    end 

    check_1 = ξᶜₐ >= 0 && ξᶜₐ <= l₀ₐ
    check_2 = ξᶜᵦ >=0 && ξᶜᵦ <= l₀ᵦ
    if check_1==1 && check_2 == 1
        dᵇˡ = norm(xᶜₐ-xᶜᵦ)
    else
        dᵇˡ = 1E9
    end
    if norm_res>tol
        dᵇˡ = 1E9
    end

    return dᵇˡ, ξᶜₐ, ξᶜᵦ

end 

function beam2beam_compute_Kc_Tc(beaminfoₐ, beaminfoᵦ, ξᶜₐ, ξᶜᵦ, pₙ, p′ₙ)

    u₁ₐ, u₂ₐ, R₁ₐ, R₂ₐ, initₐ = beaminfoₐ
    u₁ᵦ, u₂ᵦ, R₁ᵦ, R₂ᵦ, initᵦ = beaminfoᵦ

    xᶜₐ, dxᶜₐ, ddxᶜₐ  = beam2beam_get_derivatives(ξᶜₐ, u₁ₐ, u₂ₐ, R₁ₐ, R₂ₐ, initₐ)
    xᶜᵦ, dxᶜᵦ, ddxᶜᵦ  = beam2beam_get_derivatives(ξᶜᵦ, u₁ᵦ, u₂ᵦ, R₁ᵦ, R₂ᵦ, initᵦ)

    RₑH₁¹Rₑᵀₐ, RₑH₁²Rₑᵀₐ, RₑH₁³Rₑᵀₐ, RₑH₁⁴Rₑᵀₐ, RₑdH₁¹Rₑᵀₐ, RₑdH₁²Rₑᵀₐ, RₑdH₁³Rₑᵀₐ, RₑdH₁⁴Rₑᵀₐ = beam2beam_get_RₑH₁Rₑᵀ_and_RₑdH₁Rₑᵀ(ξᶜₐ, u₁ₐ, u₂ₐ, R₁ₐ, R₂ₐ, initₐ) 
    RₑH₁¹Rₑᵀᵦ, RₑH₁²Rₑᵀᵦ, RₑH₁³Rₑᵀᵦ, RₑH₁⁴Rₑᵀᵦ, RₑdH₁¹Rₑᵀᵦ, RₑdH₁²Rₑᵀᵦ, RₑdH₁³Rₑᵀᵦ, RₑdH₁⁴Rₑᵀᵦ = beam2beam_get_RₑH₁Rₑᵀ_and_RₑdH₁Rₑᵀ(ξᶜᵦ, u₁ᵦ, u₂ᵦ, R₁ᵦ, R₂ᵦ, initᵦ)

    xᶜₐᵦ = xᶜₐ-xᶜᵦ
    lₐᵦ = norm(xᶜₐᵦ)
    n = xᶜₐᵦ/norm(xᶜₐᵦ)

    kᶜ = 1E2
    Tᶜ¹ₐ = kᶜ * pₙ * RₑH₁¹Rₑᵀₐ' * n
    Tᶜ²ₐ = kᶜ * pₙ * RₑH₁²Rₑᵀₐ' * n
    Tᶜ³ₐ = kᶜ * pₙ * RₑH₁³Rₑᵀₐ' * n
    Tᶜ⁴ₐ = kᶜ * pₙ * RₑH₁⁴Rₑᵀₐ' * n

    Tᶜ¹ᵦ = kᶜ * pₙ * RₑH₁¹Rₑᵀᵦ' * n
    Tᶜ²ᵦ = kᶜ * pₙ * RₑH₁²Rₑᵀᵦ' * n
    Tᶜ³ᵦ = kᶜ * pₙ * RₑH₁³Rₑᵀᵦ' * n
    Tᶜ⁴ᵦ = kᶜ * pₙ * RₑH₁⁴Rₑᵀᵦ' * n

    Tᶜₐ = [Tᶜ¹ₐ; Tᶜ²ₐ; Tᶜ³ₐ; Tᶜ⁴ₐ]
    Tᶜᵦ = [Tᶜ¹ᵦ; Tᶜ²ᵦ; Tᶜ³ᵦ; Tᶜ⁴ᵦ]

    Nₐ = [RₑH₁¹Rₑᵀₐ; RₑH₁²Rₑᵀₐ; RₑH₁³Rₑᵀₐ; RₑH₁⁴Rₑᵀₐ]'
    Nᵦ = [RₑH₁¹Rₑᵀᵦ; RₑH₁²Rₑᵀᵦ; RₑH₁³Rₑᵀᵦ; RₑH₁⁴Rₑᵀᵦ]'
    # Nₐᵦ= [Nₐ'; -Nᵦ']'
    
    Nₐᵦ = Mat324(
        Nₐ[1,1], Nₐ[2,1], Nₐ[3,1], 
        Nₐ[1,2], Nₐ[2,2], Nₐ[3,2], 
        Nₐ[1,3], Nₐ[2,3], Nₐ[3,3], 
        Nₐ[1,4], Nₐ[2,4], Nₐ[3,4], 
        Nₐ[1,5], Nₐ[2,5], Nₐ[3,5], 
        Nₐ[1,6], Nₐ[2,6], Nₐ[3,6], 
        Nₐ[1,7], Nₐ[2,7], Nₐ[3,7], 
        Nₐ[1,8], Nₐ[2,8], Nₐ[3,8], 
        Nₐ[1,9], Nₐ[2,9], Nₐ[3,9], 
        Nₐ[1,10], Nₐ[2,10], Nₐ[3,10], 
        Nₐ[1,11], Nₐ[2,11], Nₐ[3,11], 
        Nₐ[1,12], Nₐ[2,12], Nₐ[3,12], 
        -Nᵦ[1,1], -Nᵦ[2,1], -Nᵦ[3,1], 
        -Nᵦ[1,2], -Nᵦ[2,2], -Nᵦ[3,2], 
        -Nᵦ[1,3], -Nᵦ[2,3], -Nᵦ[3,3], 
        -Nᵦ[1,4], -Nᵦ[2,4], -Nᵦ[3,4], 
        -Nᵦ[1,5],- Nᵦ[2,5], -Nᵦ[3,5], 
        -Nᵦ[1,6], -Nᵦ[2,6], -Nᵦ[3,6], 
        -Nᵦ[1,7], -Nᵦ[2,7], -Nᵦ[3,7], 
        -Nᵦ[1,8], -Nᵦ[2,8], -Nᵦ[3,8], 
        -Nᵦ[1,9], -Nᵦ[2,9], -Nᵦ[3,9], 
        -Nᵦ[1,10], -Nᵦ[2,10], -Nᵦ[3,10], 
        -Nᵦ[1,11], -Nᵦ[2,11], -Nᵦ[3,11], 
        -Nᵦ[1,12], -Nᵦ[2,12], -Nᵦ[3,12])

    dNₐ = [RₑdH₁¹Rₑᵀₐ; RₑdH₁²Rₑᵀₐ; RₑdH₁³Rₑᵀₐ; RₑdH₁⁴Rₑᵀₐ]'
    dNᵦ = [RₑdH₁¹Rₑᵀᵦ; RₑdH₁²Rₑᵀᵦ; RₑdH₁³Rₑᵀᵦ; RₑdH₁⁴Rₑᵀᵦ]'

    dgddₐᵦ = n'*Nₐᵦ
    dndₐᵦ = 1/lₐᵦ*(ID3.-n'*n)*Nₐᵦ
    dTᶜₐdₐᵦ = p′ₙ*(Nₐ'*n*dgddₐᵦ + pₙ*Nₐ'*dndₐᵦ)
    dTᶜᵦdₐᵦ = p′ₙ*(Nᵦ'*n*dgddₐᵦ + pₙ*Nᵦ'*dndₐᵦ)

    dndξₐ = 1/lₐᵦ*(ID3.-n'*n)*dxᶜₐ
    dndξᵦ = -1/lₐᵦ*(ID3.-n'*n)*dxᶜᵦ

    dTᶜₐdξₐ = pₙ*(dNₐ'*n + Nₐ'*dndξₐ)
    dTᶜₐdξᵦ= pₙ*(Nₐ'*dndξᵦ)
    dTᶜᵦdξₐ = pₙ*(Nᵦ'*dndξₐ)
    dTᶜᵦdξᵦ = pₙ*(dNᵦ'*n + Nᵦ'*dndξᵦ)

    A = Mat22(dot(dxᶜₐ, dxᶜₐ)+dot((xᶜₐ-xᶜᵦ),ddxᶜₐ), dot(dxᶜₐ,dxᶜᵦ), -dot(dxᶜₐ,dxᶜᵦ), -dot(dxᶜᵦ,dxᶜᵦ) + dot((xᶜₐ-xᶜᵦ),ddxᶜᵦ))
    aux1 = (xᶜₐ-xᶜᵦ)'*dNₐ + dxᶜₐ'*Nₐ
    aux2 = -dxᶜₐ'*Nᵦ
    aux3 = dxᶜᵦ'*Nₐ
    aux4 = (xᶜₐ-xᶜᵦ)'*dNᵦ - dxᶜᵦ'*Nᵦ
    B = Mat224(aux1[1], aux3[1], aux1[2], aux3[2], aux1[3], aux3[3], aux1[4], aux3[4], aux1[5], aux3[5], aux1[6], aux3[6],aux1[7], aux3[7], aux1[8], aux3[8], aux1[9], aux3[9], aux1[10], aux3[10],  aux1[11], aux3[11], aux1[12], aux3[12], aux2[1], aux4[1], aux2[2], aux4[2], aux2[3], aux4[3],  aux2[4], aux4[4], aux2[5], aux4[5], aux2[6], aux4[6],aux2[7], aux4[7], aux2[8], aux4[8], aux2[9], aux4[9], aux2[10], aux4[10], aux2[11], aux4[11], aux2[12], aux4[12])
    
    aux = A\(-B)

    dξₐdₐᵦ = Vec24(aux[1,1], aux[1,2], aux[1,3], aux[1,4], aux[1,5], aux[1,6], aux[1,7], aux[1,8], aux[1,9], aux[1,10], aux[1,11], aux[1,12], aux[1,13], aux[1,14], aux[1,15], aux[1,16], aux[1,17], aux[1,18], aux[1,19], aux[1,20],aux[1,21],aux[1,22],aux[1,23],aux[1,24])
    dξᵦdₐᵦ = Vec24(aux[2,1], aux[2,2], aux[2,3], aux[2,4], aux[2,5], aux[2,6], aux[2,7], aux[2,8], aux[2,9], aux[2,10], aux[2,11], aux[2,12], aux[2,13], aux[2,14], aux[2,15], aux[2,16], aux[2,17], aux[2,18], aux[2,19], aux[2,20],aux[2,21],aux[2,22],aux[2,23],aux[2,24])

    Kᶜₐ = dTᶜₐdₐᵦ + dTᶜₐdξₐ*dξₐdₐᵦ' + dTᶜₐdξᵦ*dξᵦdₐᵦ'
    Kᶜᵦ = dTᶜᵦdₐᵦ + dTᶜᵦdξₐ*dξₐdₐᵦ' + dTᶜᵦdξᵦ*dξᵦdₐᵦ'

    return Tᶜₐ, Tᶜᵦ, Kᶜₐ, Kᶜᵦ

end

function beams2beam_search_canditates(conf)

    @unpack nodes, beams = conf

    nbeams = length(beams)
    nnodes = length(nodes)
    dist_elements = 1E6*ones(Int, nbeams, nbeams)
    dist_elements_ini = 1E6*ones(Int, nbeams)
    dist_comp = 1E6*ones(Int, nbeams, nbeams)
    fact_safe = 1.1
    TNOneighbour = ones(Int, nbeams, nbeams)
    Xlist = zeros(Int, nnodes, 2)

    for (idx, b) in enumerate(beams)

        n₁ = b.node1
        n₂ = b.node2
        
        if Xlist[n₁, 1] == 0
            Xlist[n₁, 1] = idx
        else 
            Xlist[n₁, 2] = idx
        end 
        if Xlist[n₂, 1] == 0
            Xlist[n₂, 1] = idx
        else 
            Xlist[n₂, 2] = idx
        end 

    end 

    ind = findall(!iszero, Xlist[:,2])
 
    for i = eachindex(ind)
        e_neigh = Xlist[ind[i],:]
        beamₐ = e_neigh[1]
        beamᵦ = e_neigh[2]
        TNOneighbour[beamₐ,beamᵦ] = 0
        TNOneighbour[beamᵦ,beamₐ] = 0
    end

    for i = 1:nbeams
        TNOneighbour[i,i] = 0
    end

    for beamₐ in beams
        
        n₁ = beamₐ.node1
        n₂ = beamₐ.node2
        X₁, X₂ = nodes.X₀[n₁], nodes.X₀[n₂]
        u₁, u₂ = nodes.u[n₁], nodes.u[n₂]
        x₁ =  X₁ + u₁
        x₂ =  X₂ + u₂
        lₙb₁ = norm(x₂ - x₁)

        dist_elements_ini[beamₐ.ind] =  beamₐ.l₀

        xb₁ = 0.5*(x₁+x₂)

        for beamᵦ in beams[beamₐ.ind+1:end]
            
            n₁ = beamᵦ.node1
            n₂ = beamᵦ.node2
            X₁, X₂ = nodes.X₀[n₁], nodes.X₀[n₂]
            u₁, u₂ = nodes.u[n₁], nodes.u[n₂]
            x₁ =  X₁ + u₁
            x₂ =  X₂ + u₂
            lₙb₂ = norm(x₂ - x₁)
        
            xb₂ = 0.5*(x₁+x₂)
    
            dist_b₁_b₂= norm(xb₁-xb₂)
            dist_elements[beamₐ.ind,beamᵦ.ind] = dist_b₁_b₂
            dist_elements[beamᵦ.ind,beamₐ.ind] = dist_b₁_b₂
            
            dist_comp[beamₐ.ind,beamᵦ.ind] = 0.5*fact_safe*(lₙb₁+lₙb₂)
            dist_comp[beamᵦ.ind,beamₐ.ind] = 0.5*fact_safe*(lₙb₁+lₙb₂)

        end 

    end 

    check_2 = dist_elements .< dist_comp
    check_final = TNOneighbour .* check_2

    aux = findall(x->x==1, check_final) 

    candidate_elements = Vec2[]
    for i in eachindex(aux)
        push!(candidate_elements, sort([aux[i][1],aux[i][2]]))
    end 

    candidate_elements = unique(candidate_elements, dims=1)

    return candidate_elements

end 

function beam2beam_compute!(candidate_elements, conf, matrices)

    @unpack nodes, beams = conf

    for idx = eachindex(candidate_elements)

        beamₐ = beams[candidate_elements[idx][1]]
        beamᵦ = beams[candidate_elements[idx][2]]

        beaminfoₐ = get_beam_info(beamₐ, nodes)
        beaminfoᵦ = get_beam_info(beamᵦ, nodes)

        @timeit_debug "Compute distance candidates" dᵇˡ, ξᶜₐ, ξᶜᵦ =  beam2beam_compute_dist_candidate(beaminfoₐ, beaminfoᵦ)

        if dᵇˡ < 1e9

            pₙ, p′ₙ, _ = regularize_gₙ(dᵇˡ - 2*beamₐ.radius, beamₐ.radius) 

            if pₙ !=0 
                @timeit_debug "Compute force and matrix"  Tᶜₐ, Tᶜᵦ, Kᶜₐ, Kᶜᵦ =  beam2beam_compute_Kc_Tc(beaminfoₐ, beaminfoᵦ, ξᶜₐ, ξᶜᵦ, pₙ, p′ₙ)

                @timeit_debug "Assemble" begin
                    idofₐ₁ = nodes.idof_6[beamₐ.node1]
                    idofₐ₂ = nodes.idof_6[beamₐ.node2]
                    idofₐ = vcat(idofₐ₁, idofₐ₂)

                    idofᵦ₁ = nodes.idof_6[beamᵦ.node1]
                    idofᵦ₂ = nodes.idof_6[beamᵦ.node2]
                    idofᵦ = vcat(idofᵦ₁, idofᵦ₂)

                    matrices.Tᶜ[idofₐ] .+= Tᶜₐ
                    matrices.Tᶜ[idofᵦ] .+= -Tᶜᵦ

                    idofₐᵦ = vcat(idofₐ, idofᵦ)
                    matrices.K[idofₐ, idofₐᵦ] .+= -Kᶜₐ
                    matrices.K[idofᵦ, idofₐᵦ] .+= Kᶜᵦ
                end
            end 
        end 
    end
end 




