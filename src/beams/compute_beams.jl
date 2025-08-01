function compute_beams(u₁, u₂, R₁, R₂, ΔR₁, ΔR₂, u̇₁, u̇₂, ẇ₁, ẇ₂, ü₁, ü₂, ẅ₁, ẅ₂, constants, exact=true, isdynamic=true) 

    # Superscript ¹ means matrix or vector associated to u₁
    # Superscript ² means matrix or vector associated to Θ₁
    # Superscript ³ means matrix or vector associated to u₂
    # Superscript ⁴ means matrix or vector associated to Θ₂

    init, gauss_params, beamproperties = constants
    X₁, X₂, l₀, Rₑ⁰ = init
    nᴳ, ωᴳ, zᴳ = gauss_params
    @unpack K̄ⁱⁿᵗ, Jᵨ, Aᵨ, damping = beamproperties

    x₁ =  X₁ + u₁
    x₂ =  X₂ + u₂
    
    lₙ = norm(x₂ - x₁)

    ū = lₙ - l₀

    Rₑ, r¹, r³, η, Gᵀ¹, Gᵀ², Gᵀ³, Gᵀ⁴, D₃ = local_Rₑ_and_aux(x₁, x₂, R₁, R₂, Rₑ⁰[:,2], lₙ)

    R̅₁ = Rₑ' * R₁ * Rₑ⁰
    R̅₂ = Rₑ' * R₂ * Rₑ⁰

    Θ̅₁ = toangle(R̅₁)
    Θ̅₂ = toangle(R̅₂)

    if exact
        Tₛ⁻¹Θ̅₁ = Tₛ⁻¹(Θ̅₁)
        Tₛ⁻¹Θ̅₂ = Tₛ⁻¹(Θ̅₂)
    end

    P¹¹ = -Gᵀ¹ 
    P²¹ = P¹¹
    P¹² = ID3-Gᵀ²
    P²² = -Gᵀ²
    P¹⁴ = -Gᵀ⁴
    P²⁴ = ID3-Gᵀ⁴

    # B̄⁺ = [r; PEᵀ]
    B̄⁺¹ = r¹
    B̄⁺¹¹ = P¹¹ * Rₑ'
    B̄⁺²¹ = B̄⁺¹¹
    B̄⁺¹² = P¹² * Rₑ'
    B̄⁺²² = P²² * Rₑ'
    B̄⁺¹⁴ = P¹⁴ * Rₑ'
    B̄⁺²⁴ = P²⁴ * Rₑ'

    # B = B̄B̄⁺
    B¹ = B̄⁺¹
    B¹¹ = exact ?  Tₛ⁻¹Θ̅₁ * B̄⁺¹¹ : B̄⁺¹¹
    B¹² = exact ?  Tₛ⁻¹Θ̅₁ * B̄⁺¹² : B̄⁺¹²
    B¹⁴ = exact ?  Tₛ⁻¹Θ̅₁ * B̄⁺¹⁴ : B̄⁺¹⁴
    B²¹ = exact ?  Tₛ⁻¹Θ̅₂ * B̄⁺²¹ : B̄⁺²¹
    B²² = exact ?  Tₛ⁻¹Θ̅₂ * B̄⁺²² : B̄⁺²²
    B²⁴ = exact ?  Tₛ⁻¹Θ̅₂ * B̄⁺²⁴ : B̄⁺²⁴

    K̄ⁱⁿᵗū, K̄ⁱⁿᵗΘ̅, K̄ⁱⁿᵗΘ̅Θ̅ = K̄ⁱⁿᵗ

    # T̄ⁱⁿᵗ = K̄ⁱⁿᵗ D̄
    T̄ⁱⁿᵗū  = K̄ⁱⁿᵗū  * ū
    T̄ⁱⁿᵗΘ̅₁ = K̄ⁱⁿᵗΘ̅  * Θ̅₁ + K̄ⁱⁿᵗΘ̅Θ̅ * Θ̅₂
    T̄ⁱⁿᵗΘ̅₂ = K̄ⁱⁿᵗΘ̅Θ̅ * Θ̅₁ + K̄ⁱⁿᵗΘ̅  * Θ̅₂

    strain_energy = (ū*T̄ⁱⁿᵗū + dot(Θ̅₁, T̄ⁱⁿᵗΘ̅₁) + dot(Θ̅₂, T̄ⁱⁿᵗΘ̅₂))/2

    # Tⁱⁿᵗ = Bᵀ T̄ⁱⁿᵗ
    Tⁱⁿᵗ¹ = B¹'*T̄ⁱⁿᵗū + B¹¹'*T̄ⁱⁿᵗΘ̅₁ + B²¹'*T̄ⁱⁿᵗΘ̅₂
    Tⁱⁿᵗ² =             B¹²'*T̄ⁱⁿᵗΘ̅₁ + B²²'*T̄ⁱⁿᵗΘ̅₂
    Tⁱⁿᵗ³ = -Tⁱⁿᵗ¹
    Tⁱⁿᵗ⁴ =             B¹⁴'*T̄ⁱⁿᵗΘ̅₁ + B²⁴'*T̄ⁱⁿᵗΘ̅₂

    # Force
    Tⁱⁿᵗ = [Tⁱⁿᵗ¹; Tⁱⁿᵗ²; Tⁱⁿᵗ³; Tⁱⁿᵗ⁴]

    # [N̄ M̄⁺₁ M̄⁺₂] = B̄ᵀ T̄ⁱⁿᵗ
    N̄   = T̄ⁱⁿᵗū
    M̄⁺₁ = exact ? Tₛ⁻¹Θ̅₁' * T̄ⁱⁿᵗΘ̅₁  : T̄ⁱⁿᵗΘ̅₁
    M̄⁺₂ = exact ? Tₛ⁻¹Θ̅₂' * T̄ⁱⁿᵗΘ̅₂  : T̄ⁱⁿᵗΘ̅₂

    # Qₛ = Pᵀ [M̄⁺₁ M̄⁺₂]
    Qₛ¹ = P¹¹' * M̄⁺₁ + P²¹' * M̄⁺₂
    Qₛ² = P¹²' * M̄⁺₁ + P²²' * M̄⁺₂
    Qₛ⁴ = P¹⁴' * M̄⁺₁ + P²⁴' * M̄⁺₂
    
    # Q = S(Qₛ)
    Q¹ = skew(Qₛ¹)
    Q² = skew(Qₛ²)
    Q⁴ = skew(Qₛ⁴)

    a = @SVector [0, η*(M̄⁺₁[1] + M̄⁺₂[1])/lₙ + (M̄⁺₁[2] + M̄⁺₂[2])/lₙ, (M̄⁺₁[3] + M̄⁺₂[3])/lₙ]

    # DN̄ (DN̄¹¹ = DN̄³³ = -DN̄¹³ = -DN̄³¹)
    DN̄¹¹ = D₃*N̄

    #QGᵀ
    QGᵀ¹¹ = Q¹*Gᵀ¹
    QGᵀ¹² = Q¹*Gᵀ²
    QGᵀ¹⁴ = Q¹*Gᵀ⁴
    QGᵀ²² = Q²*Gᵀ²
    QGᵀ²⁴ = Q²*Gᵀ⁴
    QGᵀ⁴⁴ = Q⁴*Gᵀ⁴

    # EGa (diagonal)
    # Note: Rₑ*Ga = 0 for Θ indices because Rₑ*GᵀΘ' has only non-zero values in the first column and a = [0 ...]
    EGa¹ = Rₑ*Gᵀ¹'*a

    # EGar (EGar¹¹ = EGar³³ = -EGar³¹ = -EGar¹³)
    EGar¹¹ = EGa¹*r¹

    # Kₘ = DN̄ - EQGᵀEᵀ + EGar
    Kₘ¹¹ = DN̄¹¹ - Rₑ*QGᵀ¹¹*Rₑ' + EGar¹¹
    Kₘ¹² =      - Rₑ*QGᵀ¹²*Rₑ'
    Kₘ¹⁴ =      - Rₑ*QGᵀ¹⁴*Rₑ'
    Kₘ²² =      - Rₑ*QGᵀ²²*Rₑ'
    Kₘ²⁴ =      - Rₑ*QGᵀ²⁴*Rₑ'
    Kₘ⁴⁴ =      - Rₑ*QGᵀ⁴⁴*Rₑ'

    # K̃
    if exact
        η₁, μ₁ = compute_η_μ(Θ̅₁)
        η₂, μ₂ = compute_η_μ(Θ̅₂)

        M̄₁ = T̄ⁱⁿᵗΘ̅₁
        M̄₂ = T̄ⁱⁿᵗΘ̅₂

        K̄ₕ₁ = compute_K̄ₕ(Θ̅₁, M̄₁, Tₛ⁻¹Θ̅₁, η₁, μ₁)
        K̄ₕ₂ = compute_K̄ₕ(Θ̅₂, M̄₂, Tₛ⁻¹Θ̅₂, η₂, μ₂)
    end

    K̃¹¹ = exact ?  B̄⁺¹¹' * K̄ₕ₁ * B̄⁺¹¹  +  B̄⁺²¹' * K̄ₕ₂ * B̄⁺²¹  +  Kₘ¹¹   :   Kₘ¹¹ 
    K̃¹² = exact ?  B̄⁺¹¹' * K̄ₕ₁ * B̄⁺¹²  +  B̄⁺²¹' * K̄ₕ₂ * B̄⁺²²  +  Kₘ¹²   :   Kₘ¹² 
    K̃¹⁴ = exact ?  B̄⁺¹¹' * K̄ₕ₁ * B̄⁺¹⁴  +  B̄⁺²¹' * K̄ₕ₂ * B̄⁺²⁴  +  Kₘ¹⁴   :   Kₘ¹⁴ 
    K̃²² = exact ?  B̄⁺¹²' * K̄ₕ₁ * B̄⁺¹²  +  B̄⁺²²' * K̄ₕ₂ * B̄⁺²²  +  Kₘ²²   :   Kₘ²² 
    K̃²⁴ = exact ?  B̄⁺¹²' * K̄ₕ₁ * B̄⁺¹⁴  +  B̄⁺²²' * K̄ₕ₂ * B̄⁺²⁴  +  Kₘ²⁴   :   Kₘ²⁴ 
    K̃⁴⁴ = exact ?  B̄⁺¹⁴' * K̄ₕ₁ * B̄⁺¹⁴  +  B̄⁺²⁴' * K̄ₕ₂ * B̄⁺²⁴  +  Kₘ⁴⁴   :   Kₘ⁴⁴ 

    # BᵀKⁱⁿᵗB
    BᵀKⁱⁿᵗB¹¹ = B¹' * K̄ⁱⁿᵗū * B¹      +      B¹¹' * (K̄ⁱⁿᵗΘ̅ * B¹¹  +  K̄ⁱⁿᵗΘ̅Θ̅ * B²¹)  +  B²¹' * (K̄ⁱⁿᵗΘ̅Θ̅ * B¹¹  +  K̄ⁱⁿᵗΘ̅ * B²¹)    
    BᵀKⁱⁿᵗB¹² =                              B¹¹' * (K̄ⁱⁿᵗΘ̅ * B¹²  +  K̄ⁱⁿᵗΘ̅Θ̅ * B²²)  +  B²¹' * (K̄ⁱⁿᵗΘ̅Θ̅ * B¹²  +  K̄ⁱⁿᵗΘ̅ * B²²)     
    BᵀKⁱⁿᵗB¹⁴ =                              B¹¹' * (K̄ⁱⁿᵗΘ̅ * B¹⁴  +  K̄ⁱⁿᵗΘ̅Θ̅ * B²⁴)  +  B²¹' * (K̄ⁱⁿᵗΘ̅Θ̅ * B¹⁴  +  K̄ⁱⁿᵗΘ̅ * B²⁴)    
        
    BᵀKⁱⁿᵗB²² =                              B¹²' * (K̄ⁱⁿᵗΘ̅ * B¹²  +  K̄ⁱⁿᵗΘ̅Θ̅ * B²²)  +  B²²' * (K̄ⁱⁿᵗΘ̅Θ̅ * B¹²  +  K̄ⁱⁿᵗΘ̅ * B²²)      
    BᵀKⁱⁿᵗB²⁴ =                              B¹²' * (K̄ⁱⁿᵗΘ̅ * B¹⁴  +  K̄ⁱⁿᵗΘ̅Θ̅ * B²⁴)  +  B²²' * (K̄ⁱⁿᵗΘ̅Θ̅ * B¹⁴  +  K̄ⁱⁿᵗΘ̅ * B²⁴)      
        
    BᵀKⁱⁿᵗB⁴⁴ =                              B¹⁴' * (K̄ⁱⁿᵗΘ̅ * B¹⁴  +  K̄ⁱⁿᵗΘ̅Θ̅ * B²⁴)  +  B²⁴' * (K̄ⁱⁿᵗΘ̅Θ̅ * B¹⁴  +  K̄ⁱⁿᵗΘ̅ * B²⁴)    

    Kⁱⁿᵗ¹¹ =  BᵀKⁱⁿᵗB¹¹   +    K̃¹¹
    Kⁱⁿᵗ¹² =  BᵀKⁱⁿᵗB¹²   +    K̃¹²
    Kⁱⁿᵗ¹³ =  -Kⁱⁿᵗ¹¹
    Kⁱⁿᵗ¹⁴ =  BᵀKⁱⁿᵗB¹⁴   +    K̃¹⁴
    Kⁱⁿᵗ²² =  BᵀKⁱⁿᵗB²²   +    K̃²²
    Kⁱⁿᵗ²³ =  -Kⁱⁿᵗ¹²'
    Kⁱⁿᵗ²⁴ =  BᵀKⁱⁿᵗB²⁴   +    K̃²⁴
    Kⁱⁿᵗ³³ =  Kⁱⁿᵗ¹¹
    Kⁱⁿᵗ³⁴ =  -Kⁱⁿᵗ¹⁴
    Kⁱⁿᵗ⁴⁴ =  BᵀKⁱⁿᵗB⁴⁴   +    K̃⁴⁴

    Kⁱⁿᵗ = hcat(vcat(Kⁱⁿᵗ¹¹, Kⁱⁿᵗ¹²', Kⁱⁿᵗ¹³', Kⁱⁿᵗ¹⁴'), vcat(Kⁱⁿᵗ¹², Kⁱⁿᵗ²², Kⁱⁿᵗ²³', Kⁱⁿᵗ²⁴'), vcat(Kⁱⁿᵗ¹³, Kⁱⁿᵗ²³, Kⁱⁿᵗ³³, Kⁱⁿᵗ³⁴'), vcat(Kⁱⁿᵗ¹⁴, Kⁱⁿᵗ²⁴, Kⁱⁿᵗ³⁴, Kⁱⁿᵗ⁴⁴))

    kinetic_energy = zero(Float64)

    Tᵏ¹ = zeros(Vec3{Float64})
    Tᵏ² = zeros(Vec3{Float64})
    Tᵏ³ = zeros(Vec3{Float64})
    Tᵏ⁴ = zeros(Vec3{Float64})

    M¹¹ = zeros(Mat33{Float64})
    M¹² = zeros(Mat33{Float64})
    M¹³ = zeros(Mat33{Float64})
    M¹⁴ = zeros(Mat33{Float64})
    M²¹ = zeros(Mat33{Float64})
    M²² = zeros(Mat33{Float64})
    M²³ = zeros(Mat33{Float64})
    M²⁴ = zeros(Mat33{Float64})
    M³¹ = zeros(Mat33{Float64})
    M³² = zeros(Mat33{Float64})
    M³³ = zeros(Mat33{Float64})
    M³⁴ = zeros(Mat33{Float64})
    M⁴¹ = zeros(Mat33{Float64})
    M⁴² = zeros(Mat33{Float64})
    M⁴³ = zeros(Mat33{Float64})
    M⁴⁴ = zeros(Mat33{Float64})

    Cᵏ¹¹ = zeros(Mat33{Float64})
    Cᵏ¹² = zeros(Mat33{Float64})
    Cᵏ¹³ = zeros(Mat33{Float64})
    Cᵏ¹⁴ = zeros(Mat33{Float64})
    Cᵏ²¹ = zeros(Mat33{Float64})
    Cᵏ²² = zeros(Mat33{Float64})
    Cᵏ²³ = zeros(Mat33{Float64})
    Cᵏ²⁴ = zeros(Mat33{Float64})
    Cᵏ³¹ = zeros(Mat33{Float64})
    Cᵏ³² = zeros(Mat33{Float64})
    Cᵏ³³ = zeros(Mat33{Float64})
    Cᵏ³⁴ = zeros(Mat33{Float64})
    Cᵏ⁴¹ = zeros(Mat33{Float64})
    Cᵏ⁴² = zeros(Mat33{Float64})
    Cᵏ⁴³ = zeros(Mat33{Float64})
    Cᵏ⁴⁴ = zeros(Mat33{Float64})

    if isdynamic
        
        U̇₁ = Rₑ' * u̇₁
        U̇₂ = Rₑ' * u̇₂
        Ẇ₁ = Rₑ' * ẇ₁
        Ẇ₂ = Rₑ' * ẇ₂
        
        Ü₁ = Rₑ' * ü₁
        Ü₂ = Rₑ' * ü₂
        Ẅ₁ = Rₑ' * ẅ₁
        Ẅ₂ = Rₑ' * ẅ₂

        SU̇₁ = skew(U̇₁)
        SU̇₂ = skew(U̇₂)
        SẆ₁ = skew(Ẇ₁)
        SẆ₂ = skew(Ẇ₂)

        Ẇᵉ = Gᵀ¹ * U̇₁ + Gᵀ² * Ẇ₁ + Gᵀ³ * U̇₂ + Gᵀ⁴ * Ẇ₂
        SẆᵉ = skew(Ẇᵉ)

        rḋ = dot(r¹, u̇₁) + dot(r³, u̇₂)

        Gᵀ¹Rₑᵀ = Gᵀ¹ * Rₑ'
        Gᵀ²Rₑᵀ = Gᵀ² * Rₑ'
        Gᵀ⁴Rₑᵀ = Gᵀ⁴ * Rₑ'

        RₑG¹ = Rₑ * Gᵀ¹'
        RₑG² = Rₑ * Gᵀ²'
        RₑG⁴ = Rₑ * Gᵀ⁴'

        # cycle among the Gauss positions
        for iᴳ in 1:nᴳ

            zᴳ_iᴳ = zᴳ[iᴳ]
            ωᴳ_iᴳ = ωᴳ[iᴳ]

            ξ = l₀*(zᴳ_iᴳ+1)/2

            # Shape functions
            N₁ = 1-ξ/l₀
            N₂ = 1-N₁
            N₃ = ξ*(1-ξ/l₀)^2
            N₄ = -(1-ξ/l₀)*((ξ^2)/l₀)
            N₅ = (1-3*ξ/l₀)*(1-ξ/l₀)
            N₆ = (3*ξ/l₀-2)*(ξ/l₀)
            N₇ = N₃+N₄
            N₈ = N₅+N₆-1

            uᵗ = @SVector [0, N₃*Θ̅₁[3] + N₄*Θ̅₂[3], -N₃*Θ̅₁[2] + -N₄*Θ̅₂[2]]
            Θ̄  = @SVector [N₁*Θ̅₁[1] + N₂*Θ̅₂[1], N₅*Θ̅₁[2] + N₆*Θ̅₂[2], N₅*Θ̅₁[3] + N₆*Θ̅₂[3]]

            Suᵗ = skew(uᵗ)
            SΘ̄ = skew(Θ̄)

            R̄ = ID3 + SΘ̄

            Īᵨ = R̄*Jᵨ*R̄'

            N₇lₙ = N₇/lₙ
            N₇lₙ² = N₇lₙ/lₙ
            N₈lₙ = N₈/lₙ
            N₈lₙ² = N₈lₙ/lₙ

            P₁P¹ = @SMatrix [0 0 0; 0 N₇lₙ 0;0 0 N₇lₙ]
            P₁P² = @SMatrix [0 0 0; 0 0 N₃;0 -N₃ 0]
            P₁P³ = -P₁P¹
            P₁P⁴ = @SMatrix [0 0 0; 0 0 N₄;0 -N₄ 0]

            H₁¹ = N₁*ID3 + P₁P¹ - Suᵗ*Gᵀ¹
            H₁² =          P₁P² - Suᵗ*Gᵀ²
            H₁³ = N₂*ID3 + P₁P³ - Suᵗ*Gᵀ³
            H₁⁴ =          P₁P⁴ - Suᵗ*Gᵀ⁴

            H₂¹ = @SMatrix [0 0 0; 0  0 -N₈lₙ;0 N₈lₙ 0]
            H₂² = Diagonal(@SVector [N₁, N₅, N₅])
            H₂³ = -H₂¹
            H₂⁴ = Diagonal(@SVector [N₂, N₆, N₆])

            u̇ᵗ =  P₁P¹ * U̇₁ +  P₁P² * Ẇ₁ + P₁P³ * U̇₂ + P₁P⁴ * Ẇ₂

            Su̇ᵗ = skew(u̇ᵗ)
            
            N₇rḋ = N₇lₙ² * rḋ
            Ḣ₁¹ = Diagonal(@SVector [0, -N₇rḋ, -N₇rḋ]) - Su̇ᵗ * Gᵀ¹
            Ḣ₁² =                                      - Su̇ᵗ * Gᵀ²
            Ḣ₁⁴ =                                      - Su̇ᵗ * Gᵀ⁴

            N₈rḋ = N₈lₙ² * rḋ
            Ḣ₂¹ = @SMatrix [0 0 0; 0 0 N₈rḋ; 0 -N₈rḋ 0]

            h₁ = H₁¹ * U̇₁ + H₁² * Ẇ₁ + H₁³ * U̇₂ + H₁⁴ * Ẇ₂
            h₂ = H₂¹ * U̇₁ + H₂² * Ẇ₁ + H₂³ * U̇₂ + H₂⁴ * Ẇ₂
            Sh₁ = skew(h₁)
            Sh₂ = skew(h₂)

            C₁¹ = SẆᵉ * H₁¹ + Ḣ₁¹ - H₁¹ * SẆᵉ
            C₁² = SẆᵉ * H₁² + Ḣ₁² - H₁² * SẆᵉ
            C₁³ = -C₁¹
            C₁⁴ = SẆᵉ * H₁⁴ + Ḣ₁⁴ - H₁⁴ * SẆᵉ

            C₂¹ = SẆᵉ * H₂¹ + Ḣ₂¹ - H₂¹ * SẆᵉ
            C₂² = SẆᵉ * H₂²       - H₂² * SẆᵉ
            C₂³ = -C₂¹
            C₂⁴ = SẆᵉ * H₂⁴       - H₂⁴ * SẆᵉ

            H₁F₁ = H₁¹ * SU̇₁ + H₁² * SẆ₁ + H₁³ * SU̇₂ + H₁⁴ * SẆ₂
            H₂F₁ = H₂¹ * SU̇₁ + H₂² * SẆ₁ + H₂³ * SU̇₂ + H₂⁴ * SẆ₂

            A₁ḊrE¹ = @SMatrix [0 0 0; U̇₁[2]-U̇₂[2] 0 0; U̇₁[3]-U̇₂[3] 0 0]
            C₃¹ = -Sh₁*Gᵀ¹ + N₇lₙ²*A₁ḊrE¹ + SẆᵉ*P₁P¹ + H₁F₁*Gᵀ¹
            C₃² = -Sh₁*Gᵀ² +                  SẆᵉ*P₁P² + H₁F₁*Gᵀ²
            C₃⁴ = -Sh₁*Gᵀ⁴ +                  SẆᵉ*P₁P⁴ + H₁F₁*Gᵀ⁴

            A₂ḊrE¹ = @SMatrix [0 0 0; -U̇₁[3]+U̇₂[3] 0 0; U̇₁[2]-U̇₂[2] 0 0]
            C₄¹ = -Sh₂*Gᵀ¹ + N₈lₙ²*A₂ḊrE¹ + H₂F₁*Gᵀ¹
            C₄² = -Sh₂*Gᵀ²                  + H₂F₁*Gᵀ²
            C₄⁴ = -Sh₂*Gᵀ⁴                  + H₂F₁*Gᵀ⁴

            u̇₀ = Rₑ * h₁

            H₁Eᵀd̈ = H₁¹ * Ü₁ + H₁² * Ẅ₁ + H₁³ * Ü₂ + H₁⁴ * Ẅ₂
            C₁Eᵀḋ = C₁¹ * U̇₁ + C₁² * Ẇ₁ + C₁³ * U̇₂ + C₁⁴ * Ẇ₂
            Rₑᵀü₀ = H₁Eᵀd̈ + C₁Eᵀḋ

            Ẇ₀ = h₂
            ẇ₀ = Rₑ * Ẇ₀

            H₂Eᵀd̈ = H₂¹ * Ü₁ + H₂² * Ẅ₁ + H₂³ * Ü₂ + H₂⁴ * Ẅ₂
            C₂Eᵀḋ = C₂¹ * U̇₁ + C₂² * Ẇ₁ + C₂³ * U̇₂ + C₂⁴ * Ẇ₂
            Rₑᵀẅ₀ = H₂Eᵀd̈ + C₂Eᵀḋ

            SẆ₀ = skew(Ẇ₀)
            ĪᵨRₑᵀẅ₀ = Īᵨ*Rₑᵀẅ₀
            SẆ₀Īᵨ = SẆ₀*Īᵨ
            SẆ₀ĪᵨẆ₀ = SẆ₀Īᵨ*Ẇ₀
            ĪᵨRₑᵀẅ₀SẆ₀ĪᵨẆ₀ = ĪᵨRₑᵀẅ₀ + SẆ₀ĪᵨẆ₀
            AᵨH₁¹ᵀ = Aᵨ*H₁¹'
            AᵨH₁²ᵀ = Aᵨ*H₁²'
            AᵨH₁³ᵀ = Aᵨ*H₁³'
            AᵨH₁⁴ᵀ = Aᵨ*H₁⁴'

            Tᵏ¹G = ωᴳ_iᴳ * (AᵨH₁¹ᵀ*Rₑᵀü₀ + H₂¹'*ĪᵨRₑᵀẅ₀SẆ₀ĪᵨẆ₀)
            Tᵏ¹ += Tᵏ¹G
            Tᵏ² += ωᴳ_iᴳ * (AᵨH₁²ᵀ*Rₑᵀü₀ + H₂²'*ĪᵨRₑᵀẅ₀SẆ₀ĪᵨẆ₀)
            Tᵏ³ += -Tᵏ¹G + ωᴳ_iᴳ * Aᵨ * Rₑᵀü₀
            Tᵏ⁴ += ωᴳ_iᴳ * (AᵨH₁⁴ᵀ*Rₑᵀü₀ + H₂⁴'*ĪᵨRₑᵀẅ₀SẆ₀ĪᵨẆ₀)

            if damping>0
                Tᵈ¹G = ωᴳ_iᴳ * damping*(AᵨH₁¹ᵀ*h₁ + H₂¹'*Īᵨ*h₂)
                Tᵏ¹ += Tᵈ¹G
                Tᵏ² += ωᴳ_iᴳ * damping*(AᵨH₁²ᵀ*h₁ + H₂²'*Īᵨ*h₂)
                Tᵏ³ += -Tᵈ¹G + ωᴳ_iᴳ * Aᵨ * damping * h₁
                Tᵏ⁴ += ωᴳ_iᴳ * damping*(AᵨH₁⁴ᵀ*h₁ + H₂⁴'*Īᵨ*h₂)
            end

            M¹¹G = ωᴳ_iᴳ * (AᵨH₁¹ᵀ*H₁¹ + H₂¹'*Īᵨ*H₂¹)
            M¹²G = ωᴳ_iᴳ * (AᵨH₁¹ᵀ*H₁² + H₂¹'*Īᵨ*H₂²)
            M¹⁴G = ωᴳ_iᴳ * (AᵨH₁¹ᵀ*H₁⁴ + H₂¹'*Īᵨ*H₂⁴)

            M²²G = ωᴳ_iᴳ * (AᵨH₁²ᵀ*H₁² + H₂²'*Īᵨ*H₂²)
            M²⁴G = ωᴳ_iᴳ * (AᵨH₁²ᵀ*H₁⁴ + H₂²'*Īᵨ*H₂⁴)
            
            M⁴⁴G = ωᴳ_iᴳ * (AᵨH₁⁴ᵀ*H₁⁴ + H₂⁴'*Īᵨ*H₂⁴)

            M¹¹ += M¹¹G
            M¹² += M¹²G
            M¹³ += -M¹¹G + ωᴳ_iᴳ * AᵨH₁¹ᵀ
            M¹⁴ += M¹⁴G

            M²² += M²²G
            M²³ += -M¹²G' + ωᴳ_iᴳ * AᵨH₁²ᵀ
            M²⁴ += M²⁴G

            M³³ += M¹¹G + ωᴳ_iᴳ * Aᵨ * (ID3 - H₁¹ - H₁¹')
            M³⁴ += -M¹⁴G + ωᴳ_iᴳ * Aᵨ*H₁⁴
            
            M⁴⁴ += M⁴⁴G

            SẆ₀ĪᵨmSĪᵨẆ₀ = SẆ₀Īᵨ - skew(Īᵨ*Ẇ₀)

            AᵨC₁¹C₃¹ = Aᵨ*(C₁¹ + C₃¹)
            AᵨC₁¹C₃² = Aᵨ*(C₁² + C₃²)
            AᵨC₁¹C₃⁴ = Aᵨ*(C₁⁴ + C₃⁴)

            ĪᵨC₂¹C₄¹SẆ₀ĪᵨmSĪᵨẆ₀H₂¹ = Īᵨ*(C₂¹ + C₄¹) + SẆ₀ĪᵨmSĪᵨẆ₀ * H₂¹
            ĪᵨC₂²C₄¹SẆ₀ĪᵨmSĪᵨẆ₀H₂² = Īᵨ*(C₂² + C₄²) + SẆ₀ĪᵨmSĪᵨẆ₀ * H₂²
            ĪᵨC₂⁴C₄¹SẆ₀ĪᵨmSĪᵨẆ₀H₂⁴ = Īᵨ*(C₂⁴ + C₄⁴) + SẆ₀ĪᵨmSĪᵨẆ₀ * H₂⁴

            Cᵏ¹¹G = ωᴳ_iᴳ * ( H₁¹'*AᵨC₁¹C₃¹ + H₂¹'* ĪᵨC₂¹C₄¹SẆ₀ĪᵨmSĪᵨẆ₀H₂¹ ) 
            Cᵏ¹²G = ωᴳ_iᴳ * ( H₁¹'*AᵨC₁¹C₃² + H₂¹'* ĪᵨC₂²C₄¹SẆ₀ĪᵨmSĪᵨẆ₀H₂² ) 
            Cᵏ¹⁴G = ωᴳ_iᴳ * ( H₁¹'*AᵨC₁¹C₃⁴ + H₂¹'* ĪᵨC₂⁴C₄¹SẆ₀ĪᵨmSĪᵨẆ₀H₂⁴ ) 
 
            Cᵏ²¹G = ωᴳ_iᴳ * ( H₁²'*AᵨC₁¹C₃¹ + H₂²'* ĪᵨC₂¹C₄¹SẆ₀ĪᵨmSĪᵨẆ₀H₂¹ ) 
            Cᵏ²²G = ωᴳ_iᴳ * ( H₁²'*AᵨC₁¹C₃² + H₂²'* ĪᵨC₂²C₄¹SẆ₀ĪᵨmSĪᵨẆ₀H₂² ) 
            Cᵏ²⁴G = ωᴳ_iᴳ * ( H₁²'*AᵨC₁¹C₃⁴ + H₂²'* ĪᵨC₂⁴C₄¹SẆ₀ĪᵨmSĪᵨẆ₀H₂⁴ ) 
 
            Cᵏ⁴¹G = ωᴳ_iᴳ * ( H₁⁴'*AᵨC₁¹C₃¹ + H₂⁴'* ĪᵨC₂¹C₄¹SẆ₀ĪᵨmSĪᵨẆ₀H₂¹ ) 
            Cᵏ⁴²G = ωᴳ_iᴳ * ( H₁⁴'*AᵨC₁¹C₃² + H₂⁴'* ĪᵨC₂²C₄¹SẆ₀ĪᵨmSĪᵨẆ₀H₂² ) 
            Cᵏ⁴⁴G = ωᴳ_iᴳ * ( H₁⁴'*AᵨC₁¹C₃⁴ + H₂⁴'* ĪᵨC₂⁴C₄¹SẆ₀ĪᵨmSĪᵨẆ₀H₂⁴ ) 

            Cᵏ¹¹ += Cᵏ¹¹G
            Cᵏ¹² += Cᵏ¹²G
            Cᵏ¹⁴ += Cᵏ¹⁴G

            Cᵏ²¹ += Cᵏ²¹G
            Cᵏ²² += Cᵏ²²G
            Cᵏ²⁴ += Cᵏ²⁴G

            Cᵏ³¹ += -Cᵏ¹¹G + ωᴳ_iᴳ * AᵨC₁¹C₃¹
            Cᵏ³² += -Cᵏ¹²G + ωᴳ_iᴳ * AᵨC₁¹C₃²
            Cᵏ³⁴ += -Cᵏ¹⁴G + ωᴳ_iᴳ * AᵨC₁¹C₃⁴

            Cᵏ⁴¹ += Cᵏ⁴¹G
            Cᵏ⁴² += Cᵏ⁴²G
            Cᵏ⁴⁴ += Cᵏ⁴⁴G
            
            Īᵨᵍ = Rₑ*Īᵨ*Rₑ'
            kinetic_energy += ωᴳ_iᴳ/2 * (Aᵨ*u̇₀'*u̇₀ + ẇ₀'*Īᵨᵍ*ẇ₀)
  
        end

        l₀2 = l₀/2
        l₀2Rₑ = l₀2 * Rₑ

        Tᵏ¹ = l₀2Rₑ*Tᵏ¹
        Tᵏ² = l₀2Rₑ*Tᵏ²
        Tᵏ³ = l₀2Rₑ*Tᵏ³
        Tᵏ⁴ = l₀2Rₑ*Tᵏ⁴

        M¹¹ = l₀2Rₑ * M¹¹ * Rₑ'
        M¹² = l₀2Rₑ * M¹² * Rₑ'
        M¹³ = l₀2Rₑ * M¹³ * Rₑ'
        M¹⁴ = l₀2Rₑ * M¹⁴ * Rₑ'
        M²² = l₀2Rₑ * M²² * Rₑ'
        M²³ = l₀2Rₑ * M²³ * Rₑ'
        M²⁴ = l₀2Rₑ * M²⁴ * Rₑ'
        M³³ = l₀2Rₑ * M³³ * Rₑ'
        M³⁴ = l₀2Rₑ * M³⁴ * Rₑ'
        M⁴⁴ = l₀2Rₑ * M⁴⁴ * Rₑ'

        Cᵏ¹¹ = l₀2Rₑ * Cᵏ¹¹ * Rₑ' 
        Cᵏ¹² = l₀2Rₑ * Cᵏ¹² * Rₑ' 
        Cᵏ¹³ = -Cᵏ¹¹             
        Cᵏ¹⁴ = l₀2Rₑ * Cᵏ¹⁴ * Rₑ' 
        Cᵏ²¹ = l₀2Rₑ * Cᵏ²¹ * Rₑ' 
        Cᵏ²² = l₀2Rₑ * Cᵏ²² * Rₑ' 
        Cᵏ²³ = -Cᵏ²¹            
        Cᵏ²⁴ = l₀2Rₑ * Cᵏ²⁴ * Rₑ' 
        Cᵏ³¹ = l₀2Rₑ * Cᵏ³¹ * Rₑ'
        Cᵏ³² = l₀2Rₑ * Cᵏ³² * Rₑ' 
        Cᵏ³³ = -Cᵏ³¹              
        Cᵏ³⁴ = l₀2Rₑ * Cᵏ³⁴ * Rₑ' 
        Cᵏ⁴¹ = l₀2Rₑ * Cᵏ⁴¹ * Rₑ' 
        Cᵏ⁴² = l₀2Rₑ * Cᵏ⁴² * Rₑ' 
        Cᵏ⁴³ = -Cᵏ⁴¹
        Cᵏ⁴⁴ = l₀2Rₑ * Cᵏ⁴⁴ * Rₑ' 

        kinetic_energy = l₀2*kinetic_energy

        if damping>0
            Cᵏ¹¹ += damping*M¹¹
            Cᵏ¹² += damping*M¹²
            Cᵏ¹³ += damping*M¹³
            Cᵏ¹⁴ += damping*M¹⁴
            Cᵏ²¹ += damping*M¹²'
            Cᵏ²² += damping*M²²
            Cᵏ²³ += damping*M²³
            Cᵏ²⁴ += damping*M²⁴
            Cᵏ³¹ += damping*M¹³'
            Cᵏ³² += damping*M²³'
            Cᵏ³³ += damping*M³³
            Cᵏ³⁴ += damping*M³⁴
            Cᵏ⁴¹ += damping*M¹⁴'
            Cᵏ⁴² += damping*M²⁴'
            Cᵏ⁴³ += damping*M³⁴'
            Cᵏ⁴⁴ += damping*M⁴⁴
        end

        if exact
            
            Θ₁ = toangle(ΔR₁)
            Θ₂ = toangle(ΔR₂)
            
            Tₛ⁻¹Θ₁ = Tₛ⁻¹(Θ₁)
            Tₛ⁻¹Θ₂ = Tₛ⁻¹(Θ₂)

            M²¹ = M¹²' * Tₛ⁻¹Θ₁'
            M²² = M²²  * Tₛ⁻¹Θ₁' 
            M²³ = M²³  * Tₛ⁻¹Θ₁'
            M²⁴ = M²⁴  * Tₛ⁻¹Θ₁'
            M³¹ = M¹³'
            M³² = M²³'
            M⁴¹ = M¹⁴' * Tₛ⁻¹Θ₂'
            M⁴² = M²⁴' * Tₛ⁻¹Θ₂' 
            M⁴³ = M³⁴' * Tₛ⁻¹Θ₂'
            M⁴⁴ = M⁴⁴  * Tₛ⁻¹Θ₂'

            Cᵏ²¹ = Cᵏ²¹ * Tₛ⁻¹Θ₁'
            Cᵏ²² = Cᵏ²² * Tₛ⁻¹Θ₁' 
            Cᵏ²³ = Cᵏ²³ * Tₛ⁻¹Θ₁'
            Cᵏ²⁴ = Cᵏ²⁴ * Tₛ⁻¹Θ₁'
            Cᵏ⁴¹ = Cᵏ⁴¹ * Tₛ⁻¹Θ₂'
            Cᵏ⁴² = Cᵏ⁴² * Tₛ⁻¹Θ₂' 
            Cᵏ⁴³ = Cᵏ⁴³ * Tₛ⁻¹Θ₂'
            Cᵏ⁴⁴ = Cᵏ⁴⁴ * Tₛ⁻¹Θ₂'

        else

            M²¹ = M¹²'
            M³¹ = M¹³'
            M³² = M²³'
            M⁴¹ = M¹⁴' 
            M⁴² = M²⁴'
            M⁴³ = M³⁴'

        end
        
    end

    Tᵏ = [Tᵏ¹; Tᵏ²; Tᵏ³; Tᵏ⁴]
    M = vcat(
        hcat(M¹¹, M¹², M¹³, M¹⁴),
        hcat(M²¹, M²², M²³, M²⁴),
        hcat(M³¹, M³², M³³, M³⁴),
        hcat(M⁴¹, M⁴², M⁴³, M⁴⁴)
    )
    
    Cᵏ = vcat(
        hcat(Cᵏ¹¹, Cᵏ¹², Cᵏ¹³, Cᵏ¹⁴),
        hcat(Cᵏ²¹, Cᵏ²², Cᵏ²³, Cᵏ²⁴),
        hcat(Cᵏ³¹, Cᵏ³², Cᵏ³³, Cᵏ³⁴),
        hcat(Cᵏ⁴¹, Cᵏ⁴², Cᵏ⁴³, Cᵏ⁴⁴)
    )

    return strain_energy, kinetic_energy, Tⁱⁿᵗ, Tᵏ, Kⁱⁿᵗ, M, Cᵏ

end


