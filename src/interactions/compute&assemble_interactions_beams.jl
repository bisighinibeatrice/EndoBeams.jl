function compute_contact_beams(u₁, u₂, R₁, R₂, ΔR₁, ΔR₂, u̇₁, u̇₂, ẇ₁, ẇ₂, constants, exact=true) 
    
    init, gauss_params, radius_beam, properties_contact, master_surface = constants 
    X₁, X₂, l₀, Rₑ⁰ = init 
    nᴳ, ωᴳ, zᴳ = gauss_params 
    @unpack kₙ, μ, εᵗ, ηₙ, kₜ, ηₜ = properties_contact
    
    x₁ =  X₁ + u₁
    x₂ =  X₂ + u₂
    
    lₙ = norm(x₂ - x₁)
    
    Rₑ, _, r³, _, Gᵀ¹, Gᵀ², Gᵀ³, Gᵀ⁴, _ = local_Rₑ_and_aux(x₁, x₂, R₁, R₂, Rₑ⁰[:,2], lₙ)
    
    R̅₁ = Rₑ' * R₁ * Rₑ⁰
    R̅₂ = Rₑ' * R₂ * Rₑ⁰
    
    Θ̅₁ = toangle(R̅₁)
    Θ̅₂ = toangle(R̅₂)
    
    U̇₁ = Rₑ' * u̇₁
    U̇₂ = Rₑ' * u̇₂
    Ẇ₁ = Rₑ' * ẇ₁
    Ẇ₂ = Rₑ' * ẇ₂
    
    SU̇₁ = skew(U̇₁)
    SU̇₂ = skew(U̇₂)
    SẆ₁ = skew(Ẇ₁)
    SẆ₂ = skew(Ẇ₂)
    
    Gᵀ¹Rₑᵀ = Gᵀ¹ * Rₑ'
    Gᵀ²Rₑᵀ = Gᵀ² * Rₑ'
    Gᵀ⁴Rₑᵀ = Gᵀ⁴ * Rₑ'
    
    RₑG¹ = Rₑ * Gᵀ¹'
    RₑG² = Rₑ * Gᵀ²'
    RₑG⁴ = Rₑ * Gᵀ⁴'
    
    contact_energy = zero(Float64)   
    
    Tᶜ¹ = zeros(Vec3{Float64}); Tᶜ² = zeros(Vec3{Float64}); Tᶜ³ = zeros(Vec3{Float64}); Tᶜ⁴ = zeros(Vec3{Float64})
    
    Kᶜ¹¹ = zeros(Mat33{Float64}); Kᶜ¹² = zeros(Mat33{Float64}); Kᶜ¹³ = zeros(Mat33{Float64}); Kᶜ¹⁴ = zeros(Mat33{Float64})
    Kᶜ²¹ = zeros(Mat33{Float64}); Kᶜ²² = zeros(Mat33{Float64}); Kᶜ²³ = zeros(Mat33{Float64}); Kᶜ²⁴ = zeros(Mat33{Float64}) 
    Kᶜ³¹ = zeros(Mat33{Float64}); Kᶜ³² = zeros(Mat33{Float64}); Kᶜ³³ = zeros(Mat33{Float64}); Kᶜ³⁴ = zeros(Mat33{Float64})
    Kᶜ⁴¹ = zeros(Mat33{Float64}); Kᶜ⁴² = zeros(Mat33{Float64}); Kᶜ⁴³ = zeros(Mat33{Float64}); Kᶜ⁴⁴ = zeros(Mat33{Float64})
    
    Cᶜ¹¹ = zeros(Mat33{Float64}); Cᶜ¹² = zeros(Mat33{Float64}); Cᶜ¹³ = zeros(Mat33{Float64}); Cᶜ¹⁴ = zeros(Mat33{Float64})
    Cᶜ²¹ = zeros(Mat33{Float64}); Cᶜ²² = zeros(Mat33{Float64}); Cᶜ²³ = zeros(Mat33{Float64}); Cᶜ²⁴ = zeros(Mat33{Float64})
    Cᶜ³¹ = zeros(Mat33{Float64}); Cᶜ³² = zeros(Mat33{Float64}); Cᶜ³³ = zeros(Mat33{Float64}); Cᶜ³⁴ = zeros(Mat33{Float64})
    Cᶜ⁴¹ = zeros(Mat33{Float64}); Cᶜ⁴² = zeros(Mat33{Float64}); Cᶜ⁴³ = zeros(Mat33{Float64}); Cᶜ⁴⁴ = zeros(Mat33{Float64})
    
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
        Suᵗ = skew(uᵗ)
        xᴳ = N₁*x₁ + N₂*x₂ + Rₑ*uᵗ        
        
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
        
        h₁ = H₁¹ * U̇₁ + H₁² * Ẇ₁ + H₁³ * U̇₂ + H₁⁴ * Ẇ₂
        Sh₁ = skew(h₁)
        
        H₁F₁ = H₁¹ * SU̇₁ + H₁² * SẆ₁ + H₁³ * SU̇₂ + H₁⁴ * SẆ₂
        u̇₀ = Rₑ * h₁
        
        if incontact_beams(xᴳ, master_surface, radius_beam)
            
            gₙ, ∂gₙ∂x, ∂²gₙ∂x² = contact_gap_beams(xᴳ, master_surface, radius_beam) 
            
            pₙ, p′ₙ, Πₑ = regularize_gap_penalty_beams(gₙ, radius_beam)
            ηₙ, η′ₙ = smoothstep_transition_beams(ηₙ, gₙ, radius_beam)
            
            u̇ₙ_mag = dot(u̇₀, ∂gₙ∂x) #u̇₀
            u̇ₙ = u̇ₙ_mag*∂gₙ∂x
            u̇ₜ = u̇₀ - u̇ₙ
            u̇ₜ² = dot(u̇ₜ, u̇ₜ)
            
            contact_energy += ωᴳ_iᴳ*kₙ*Πₑ
            
            𝓯ⁿ = kₙ * pₙ * ∂gₙ∂x - ηₙ * u̇ₙ
            
            μʳᵉᵍ = μ/sqrt(u̇ₜ²+εᵗ)
            𝓯ᵗ = - kₙ * pₙ * μʳᵉᵍ * u̇ₜ   
            
            𝓯ᶜ = 𝓯ⁿ + 𝓯ᵗ
            
            𝓕ᶜ = Rₑ' * 𝓯ᶜ
            
            RₑH₁¹Rₑᵀ = Rₑ * H₁¹ * Rₑ'
            RₑH₁²Rₑᵀ = Rₑ * H₁² * Rₑ'
            RₑH₁³Rₑᵀ = Rₑ * H₁³ * Rₑ'
            RₑH₁⁴Rₑᵀ = Rₑ * H₁⁴ * Rₑ'
            
            Tᶜ¹ += ωᴳ_iᴳ * (RₑH₁¹Rₑᵀ' * 𝓯ᶜ)
            Tᶜ² += ωᴳ_iᴳ * (RₑH₁²Rₑᵀ' * 𝓯ᶜ)
            Tᶜ³ += ωᴳ_iᴳ * (RₑH₁³Rₑᵀ' * 𝓯ᶜ)
            Tᶜ⁴ += ωᴳ_iᴳ * (RₑH₁⁴Rₑᵀ' * 𝓯ᶜ)
            
            ŜH₁ᵀ𝓕ᶜ¹ = skew(H₁¹' * 𝓕ᶜ)
            ŜH₁ᵀ𝓕ᶜ² = skew(H₁²' * 𝓕ᶜ)
            ŜH₁ᵀ𝓕ᶜ³ = skew(H₁³' * 𝓕ᶜ)
            ŜH₁ᵀ𝓕ᶜ⁴ = skew(H₁⁴' * 𝓕ᶜ)
            
            RₑŜH₁ᵀ𝓕ᶜ¹ = Rₑ * ŜH₁ᵀ𝓕ᶜ¹
            t₁¹¹ = -RₑŜH₁ᵀ𝓕ᶜ¹ * Gᵀ¹Rₑᵀ
            t₁¹² = -RₑŜH₁ᵀ𝓕ᶜ¹ * Gᵀ²Rₑᵀ
            t₁¹³ = -t₁¹¹
            t₁¹⁴ = -RₑŜH₁ᵀ𝓕ᶜ¹ * Gᵀ⁴Rₑᵀ
            
            RₑŜH₁ᵀ𝓕ᶜ² = Rₑ * ŜH₁ᵀ𝓕ᶜ²
            t₁²¹ = -RₑŜH₁ᵀ𝓕ᶜ² * Gᵀ¹Rₑᵀ
            t₁²² = -RₑŜH₁ᵀ𝓕ᶜ² * Gᵀ²Rₑᵀ
            t₁²³ = -t₁²¹
            t₁²⁴ = -RₑŜH₁ᵀ𝓕ᶜ² * Gᵀ⁴Rₑᵀ
            
            RₑŜH₁ᵀ𝓕ᶜ³ = Rₑ * ŜH₁ᵀ𝓕ᶜ³
            t₁³¹ = -RₑŜH₁ᵀ𝓕ᶜ³ * Gᵀ¹Rₑᵀ
            t₁³² = -RₑŜH₁ᵀ𝓕ᶜ³ * Gᵀ²Rₑᵀ
            t₁³³ = -t₁³¹
            t₁³⁴ = -RₑŜH₁ᵀ𝓕ᶜ³ * Gᵀ⁴Rₑᵀ
            
            RₑŜH₁ᵀ𝓕ᶜ⁴ = Rₑ * ŜH₁ᵀ𝓕ᶜ⁴
            t₁⁴¹ = -RₑŜH₁ᵀ𝓕ᶜ⁴ * Gᵀ¹Rₑᵀ
            t₁⁴² = -RₑŜH₁ᵀ𝓕ᶜ⁴ * Gᵀ²Rₑᵀ
            t₁⁴³ = -t₁⁴¹
            t₁⁴⁴ = -RₑŜH₁ᵀ𝓕ᶜ⁴ * Gᵀ⁴Rₑᵀ
            
            v₁ = r³
            A₁ᵀ𝓕ᶜr₁₁ = @SMatrix[0 0 0; 𝓕ᶜ[2]*v₁[1] 𝓕ᶜ[2]*v₁[2] 𝓕ᶜ[2]*v₁[3]; 𝓕ᶜ[3]*v₁[1] 𝓕ᶜ[3]*v₁[2] 𝓕ᶜ[3]*v₁[3]]
            
            S𝓕ᶜ = skew(𝓕ᶜ)
            
            S𝓕ᶜP₁P¹Rₑᵀ = S𝓕ᶜ * P₁P¹ * Rₑ' 
            S𝓕ᶜP₁P²Rₑᵀ = S𝓕ᶜ * P₁P² * Rₑ'
            S𝓕ᶜP₁P⁴Rₑᵀ = S𝓕ᶜ * P₁P⁴ * Rₑ'
            
            t₂¹¹ = N₇lₙ² * Rₑ * A₁ᵀ𝓕ᶜr₁₁ - RₑG¹ * S𝓕ᶜP₁P¹Rₑᵀ 
            t₂¹² =                       - RₑG¹ * S𝓕ᶜP₁P²Rₑᵀ
            t₂¹³ = -t₂¹¹
            t₂¹⁴ =                       - RₑG¹ * S𝓕ᶜP₁P⁴Rₑᵀ
            
            t₂²¹ =                       - RₑG² * S𝓕ᶜP₁P¹Rₑᵀ
            t₂²² =                       - RₑG² * S𝓕ᶜP₁P²Rₑᵀ
            t₂²³ = -t₂²¹
            t₂²⁴ =                       - RₑG² * S𝓕ᶜP₁P⁴Rₑᵀ
            
            t₂³¹ = t₂¹³
            t₂³² = -t₂¹²
            t₂³³ = t₂¹¹
            t₂³⁴ = -t₂¹⁴
            
            t₂⁴¹ =                       - RₑG⁴ * S𝓕ᶜP₁P¹Rₑᵀ
            t₂⁴² =                       - RₑG⁴ * S𝓕ᶜP₁P²Rₑᵀ
            t₂⁴³ = -t₂⁴¹
            t₂⁴⁴ =                       - RₑG⁴ * S𝓕ᶜP₁P⁴Rₑᵀ
            
            RₑH₁ᵀ¹S𝓕ᶜ = Rₑ * H₁¹' * S𝓕ᶜ
            RₑH₁ᵀ²S𝓕ᶜ = Rₑ * H₁²' * S𝓕ᶜ
            RₑH₁ᵀ³S𝓕ᶜ = Rₑ * H₁³' * S𝓕ᶜ
            RₑH₁ᵀ⁴S𝓕ᶜ = Rₑ * H₁⁴' * S𝓕ᶜ
            
            t₃¹¹ = RₑH₁ᵀ¹S𝓕ᶜ * Gᵀ¹Rₑᵀ
            t₃¹² = RₑH₁ᵀ¹S𝓕ᶜ * Gᵀ²Rₑᵀ
            t₃¹³ = -t₃¹¹
            t₃¹⁴ = RₑH₁ᵀ¹S𝓕ᶜ * Gᵀ⁴Rₑᵀ
            
            t₃²¹ = RₑH₁ᵀ²S𝓕ᶜ * Gᵀ¹Rₑᵀ
            t₃²² = RₑH₁ᵀ²S𝓕ᶜ * Gᵀ²Rₑᵀ
            t₃²³ = -t₃²¹
            t₃²⁴ = RₑH₁ᵀ²S𝓕ᶜ * Gᵀ⁴Rₑᵀ
            
            t₃³¹ = RₑH₁ᵀ³S𝓕ᶜ * Gᵀ¹Rₑᵀ
            t₃³² = RₑH₁ᵀ³S𝓕ᶜ * Gᵀ²Rₑᵀ
            t₃³³ = -t₃³¹
            t₃³⁴ = RₑH₁ᵀ³S𝓕ᶜ * Gᵀ⁴Rₑᵀ
            
            t₃⁴¹ = RₑH₁ᵀ⁴S𝓕ᶜ * Gᵀ¹Rₑᵀ
            t₃⁴² = RₑH₁ᵀ⁴S𝓕ᶜ * Gᵀ²Rₑᵀ
            t₃⁴³ = -t₃⁴¹
            t₃⁴⁴ = RₑH₁ᵀ⁴S𝓕ᶜ * Gᵀ⁴Rₑᵀ
            
            aux = dot(∂gₙ∂x, u̇₀) * ∂²gₙ∂x² + ∂gₙ∂x * (∂²gₙ∂x² * u̇₀)'
            
            𝓐₁¹ =  aux * RₑH₁¹Rₑᵀ
            𝓐₁² =  aux * RₑH₁²Rₑᵀ
            𝓐₁³ =  aux * RₑH₁³Rₑᵀ
            𝓐₁⁴ =  aux * RₑH₁⁴Rₑᵀ
            
            RₑSh₁ = Rₑ * Sh₁
            𝓐₂¹ = - RₑSh₁ * Gᵀ¹Rₑᵀ
            𝓐₂² = - RₑSh₁ * Gᵀ²Rₑᵀ
            𝓐₂⁴ = - RₑSh₁ * Gᵀ⁴Rₑᵀ
            
            A₁Ḋr¹ = @SMatrix [0 0 0; v₁[1]*(U̇₁[2]-U̇₂[2]) v₁[2]*(U̇₁[2]-U̇₂[2]) v₁[3]*(U̇₁[2]-U̇₂[2]);v₁[1]*(U̇₁[3]-U̇₂[3]) v₁[2]*(U̇₁[3]-U̇₂[3]) v₁[3]*(U̇₁[3]-U̇₂[3])]
            
            RₑSGᵀḊ = Rₑ * skew(Gᵀ¹ * U̇₁ + Gᵀ² * Ẇ₁ + Gᵀ³ * U̇₂ + Gᵀ⁴ * Ẇ₂)
            𝓐₃¹ = N₇lₙ² * Rₑ * A₁Ḋr¹ + RₑSGᵀḊ * P₁P¹ * Rₑ'
            𝓐₃² =                      RₑSGᵀḊ * P₁P² * Rₑ'
            𝓐₃⁴ =                      RₑSGᵀḊ * P₁P⁴ * Rₑ'
            
            RₑH₁SḊ = Rₑ*H₁F₁
            𝓐₄¹ = RₑH₁SḊ * Gᵀ¹Rₑᵀ
            𝓐₄² = RₑH₁SḊ * Gᵀ²Rₑᵀ
            𝓐₄⁴ = RₑH₁SḊ * Gᵀ⁴Rₑᵀ
            
            ∂u̇₀∂d¹ = 𝓐₂¹ + 𝓐₃¹ + 𝓐₄¹
            ∂u̇₀∂d² = 𝓐₂² + 𝓐₃² + 𝓐₄²
            ∂u̇₀∂d³ = -∂u̇₀∂d¹
            ∂u̇₀∂d⁴ = 𝓐₂⁴ + 𝓐₃⁴ + 𝓐₄⁴
            
            nn = ∂gₙ∂x*∂gₙ∂x'
            
            ∂u̇ₙ∂d¹ = 𝓐₁¹ + nn * ∂u̇₀∂d¹
            ∂u̇ₙ∂d² = 𝓐₁² + nn * ∂u̇₀∂d²
            ∂u̇ₙ∂d³ = 𝓐₁³ + nn * ∂u̇₀∂d³
            ∂u̇ₙ∂d⁴ = 𝓐₁⁴ + nn * ∂u̇₀∂d⁴
            
            # u̇ₜ = u̇₀ - u̇ₙ
            ∂u̇ₜ∂d¹ = ∂u̇₀∂d¹ - ∂u̇ₙ∂d¹
            ∂u̇ₜ∂d² = ∂u̇₀∂d² - ∂u̇ₙ∂d²
            ∂u̇ₜ∂d³ = ∂u̇₀∂d³ - ∂u̇ₙ∂d³
            ∂u̇ₜ∂d⁴ = ∂u̇₀∂d⁴ - ∂u̇ₙ∂d⁴
            
            ∂u̇ₙ∂ḋ¹ = nn * RₑH₁¹Rₑᵀ
            ∂u̇ₙ∂ḋ² = nn * RₑH₁²Rₑᵀ
            ∂u̇ₙ∂ḋ³ = nn * RₑH₁³Rₑᵀ
            ∂u̇ₙ∂ḋ⁴ = nn * RₑH₁⁴Rₑᵀ
            
            ∂u̇₀∂ḋ¹ = RₑH₁¹Rₑᵀ
            ∂u̇₀∂ḋ² = RₑH₁²Rₑᵀ
            ∂u̇₀∂ḋ³ = RₑH₁³Rₑᵀ
            ∂u̇₀∂ḋ⁴ = RₑH₁⁴Rₑᵀ
            
            # u̇ₜ = u̇₀ - u̇ₙ
            ∂u̇ₜ∂ḋ¹ = ∂u̇₀∂ḋ¹ - ∂u̇ₙ∂ḋ¹
            ∂u̇ₜ∂ḋ² = ∂u̇₀∂ḋ² - ∂u̇ₙ∂ḋ²
            ∂u̇ₜ∂ḋ³ = ∂u̇₀∂ḋ³ - ∂u̇ₙ∂ḋ³
            ∂u̇ₜ∂ḋ⁴ = ∂u̇₀∂ḋ⁴ - ∂u̇ₙ∂ḋ⁴
            
            aux = kₙ*p′ₙ*nn + kₙ*pₙ*∂²gₙ∂x² - η′ₙ*u̇ₙ*∂gₙ∂x'
            Kᶠⁿ¹ = aux * RₑH₁¹Rₑᵀ - ηₙ*∂u̇ₙ∂d¹
            Kᶠⁿ² = aux * RₑH₁²Rₑᵀ - ηₙ*∂u̇ₙ∂d²
            Kᶠⁿ³ = aux * RₑH₁³Rₑᵀ - ηₙ*∂u̇ₙ∂d³
            Kᶠⁿ⁴ = aux * RₑH₁⁴Rₑᵀ - ηₙ*∂u̇ₙ∂d⁴
            
            Cᶠⁿ¹ = - ηₙ * nn * ∂u̇₀∂ḋ¹
            Cᶠⁿ² = - ηₙ * nn * ∂u̇₀∂ḋ²
            Cᶠⁿ³ = - ηₙ * nn * ∂u̇₀∂ḋ³
            Cᶠⁿ⁴ = - ηₙ * nn * ∂u̇₀∂ḋ⁴
            
            Kᶠᵗ¹ = zeros(Mat33)
            Kᶠᵗ² = zeros(Mat33)
            Kᶠᵗ³ = zeros(Mat33)
            Kᶠᵗ⁴ = zeros(Mat33)
            
            Cᶠᵗ¹ = zeros(Mat33)
            Cᶠᵗ² = zeros(Mat33)
            Cᶠᵗ³ = zeros(Mat33)
            Cᶠᵗ⁴ = zeros(Mat33)
            
            aux1 = - μʳᵉᵍ * kₙ * p′ₙ * u̇ₜ * ∂gₙ∂x'
            aux2 = - μʳᵉᵍ * kₙ * pₙ*(ID3 - 1/(u̇ₜ²+εᵗ)*u̇ₜ*u̇ₜ')
            Kᶠᵗ¹ = aux1 * RₑH₁¹Rₑᵀ + aux2 * ∂u̇ₜ∂d¹
            Kᶠᵗ² = aux1 * RₑH₁²Rₑᵀ + aux2 * ∂u̇ₜ∂d²
            Kᶠᵗ³ = aux1 * RₑH₁³Rₑᵀ + aux2 * ∂u̇ₜ∂d³
            Kᶠᵗ⁴ = aux1 * RₑH₁⁴Rₑᵀ + aux2 * ∂u̇ₜ∂d⁴
            
            Cᶠᵗ¹ = aux2 * ∂u̇ₜ∂ḋ¹
            Cᶠᵗ² = aux2 * ∂u̇ₜ∂ḋ²
            Cᶠᵗ³ = aux2 * ∂u̇ₜ∂ḋ³
            Cᶠᵗ⁴ = aux2 * ∂u̇ₜ∂ḋ⁴
            
            Kᶠᶜ¹ = Kᶠⁿ¹ + Kᶠᵗ¹
            Kᶠᶜ² = Kᶠⁿ² + Kᶠᵗ²
            Kᶠᶜ³ = Kᶠⁿ³ + Kᶠᵗ³
            Kᶠᶜ⁴ = Kᶠⁿ⁴ + Kᶠᵗ⁴
            
            t₄¹¹ = RₑH₁¹Rₑᵀ' * Kᶠᶜ¹
            t₄¹² = RₑH₁¹Rₑᵀ' * Kᶠᶜ²
            t₄¹³ = RₑH₁¹Rₑᵀ' * Kᶠᶜ³
            t₄¹⁴ = RₑH₁¹Rₑᵀ' * Kᶠᶜ⁴
            
            t₄²¹ = RₑH₁²Rₑᵀ' * Kᶠᶜ¹
            t₄²² = RₑH₁²Rₑᵀ' * Kᶠᶜ²
            t₄²³ = RₑH₁²Rₑᵀ' * Kᶠᶜ³
            t₄²⁴ = RₑH₁²Rₑᵀ' * Kᶠᶜ⁴
            
            t₄³¹ = RₑH₁³Rₑᵀ' * Kᶠᶜ¹
            t₄³² = RₑH₁³Rₑᵀ' * Kᶠᶜ²
            t₄³³ = RₑH₁³Rₑᵀ' * Kᶠᶜ³
            t₄³⁴ = RₑH₁³Rₑᵀ' * Kᶠᶜ⁴
            
            t₄⁴¹ = RₑH₁⁴Rₑᵀ' * Kᶠᶜ¹
            t₄⁴² = RₑH₁⁴Rₑᵀ' * Kᶠᶜ²
            t₄⁴³ = RₑH₁⁴Rₑᵀ' * Kᶠᶜ³
            t₄⁴⁴ = RₑH₁⁴Rₑᵀ' * Kᶠᶜ⁴
            
            Cᶠᶜ¹ = Cᶠⁿ¹ + Cᶠᵗ¹
            Cᶠᶜ² = Cᶠⁿ² + Cᶠᵗ²
            Cᶠᶜ³ = Cᶠⁿ³ + Cᶠᵗ³
            Cᶠᶜ⁴ = Cᶠⁿ⁴ + Cᶠᵗ⁴
            
            Kᶜ¹¹ +=  ωᴳ_iᴳ * (t₁¹¹ + t₂¹¹ + t₃¹¹ + t₄¹¹)
            Kᶜ¹² +=  ωᴳ_iᴳ * (t₁¹² + t₂¹² + t₃¹² + t₄¹²)
            Kᶜ¹³ +=  ωᴳ_iᴳ * (t₁¹³ + t₂¹³ + t₃¹³ + t₄¹³)
            Kᶜ¹⁴ +=  ωᴳ_iᴳ * (t₁¹⁴ + t₂¹⁴ + t₃¹⁴ + t₄¹⁴)
            
            Kᶜ²¹ +=  ωᴳ_iᴳ * (t₁²¹ + t₂²¹ + t₃²¹ + t₄²¹)
            Kᶜ²² +=  ωᴳ_iᴳ * (t₁²² + t₂²² + t₃²² + t₄²²)
            Kᶜ²³ +=  ωᴳ_iᴳ * (t₁²³ + t₂²³ + t₃²³ + t₄²³)
            Kᶜ²⁴ +=  ωᴳ_iᴳ * (t₁²⁴ + t₂²⁴ + t₃²⁴ + t₄²⁴)
            
            Kᶜ³¹ +=  ωᴳ_iᴳ * (t₁³¹ + t₂³¹ + t₃³¹ + t₄³¹)
            Kᶜ³² +=  ωᴳ_iᴳ * (t₁³² + t₂³² + t₃³² + t₄³²)
            Kᶜ³³ +=  ωᴳ_iᴳ * (t₁³³ + t₂³³ + t₃³³ + t₄³³)
            Kᶜ³⁴ +=  ωᴳ_iᴳ * (t₁³⁴ + t₂³⁴ + t₃³⁴ + t₄³⁴)
            
            Kᶜ⁴¹ +=  ωᴳ_iᴳ * (t₁⁴¹ + t₂⁴¹ + t₃⁴¹ + t₄⁴¹)
            Kᶜ⁴² +=  ωᴳ_iᴳ * (t₁⁴² + t₂⁴² + t₃⁴² + t₄⁴²)
            Kᶜ⁴³ +=  ωᴳ_iᴳ * (t₁⁴³ + t₂⁴³ + t₃⁴³ + t₄⁴³)
            Kᶜ⁴⁴ +=  ωᴳ_iᴳ * (t₁⁴⁴ + t₂⁴⁴ + t₃⁴⁴ + t₄⁴⁴)
            
            Cᶜ¹¹ +=  ωᴳ_iᴳ * RₑH₁¹Rₑᵀ' * Cᶠᶜ¹
            Cᶜ¹² +=  ωᴳ_iᴳ * RₑH₁¹Rₑᵀ' * Cᶠᶜ²
            Cᶜ¹³ +=  ωᴳ_iᴳ * RₑH₁¹Rₑᵀ' * Cᶠᶜ³
            Cᶜ¹⁴ +=  ωᴳ_iᴳ * RₑH₁¹Rₑᵀ' * Cᶠᶜ⁴
            
            Cᶜ²¹ +=  ωᴳ_iᴳ * RₑH₁²Rₑᵀ' * Cᶠᶜ¹
            Cᶜ²² +=  ωᴳ_iᴳ * RₑH₁²Rₑᵀ' * Cᶠᶜ²
            Cᶜ²³ +=  ωᴳ_iᴳ * RₑH₁²Rₑᵀ' * Cᶠᶜ³
            Cᶜ²⁴ +=  ωᴳ_iᴳ * RₑH₁²Rₑᵀ' * Cᶠᶜ⁴
            
            Cᶜ³¹ +=  ωᴳ_iᴳ * RₑH₁³Rₑᵀ' * Cᶠᶜ¹
            Cᶜ³² +=  ωᴳ_iᴳ * RₑH₁³Rₑᵀ' * Cᶠᶜ²
            Cᶜ³³ +=  ωᴳ_iᴳ * RₑH₁³Rₑᵀ' * Cᶠᶜ³
            Cᶜ³⁴ +=  ωᴳ_iᴳ * RₑH₁³Rₑᵀ' * Cᶠᶜ⁴
            
            Cᶜ⁴¹ +=  ωᴳ_iᴳ * RₑH₁⁴Rₑᵀ' * Cᶠᶜ¹
            Cᶜ⁴² +=  ωᴳ_iᴳ * RₑH₁⁴Rₑᵀ' * Cᶠᶜ²
            Cᶜ⁴³ +=  ωᴳ_iᴳ * RₑH₁⁴Rₑᵀ' * Cᶠᶜ³
            Cᶜ⁴⁴ +=  ωᴳ_iᴳ * RₑH₁⁴Rₑᵀ' * Cᶠᶜ⁴
            
        end
        
    end
    
    l₀2 = l₀/2
    l₀2Rₑ = l₀2 * Rₑ
    
    Tᶜ¹ = l₀2*Tᶜ¹
    Tᶜ² = l₀2*Tᶜ²
    Tᶜ³ = l₀2*Tᶜ³
    Tᶜ⁴ = l₀2*Tᶜ⁴
    
    Cᶜ¹¹ = l₀2 * Cᶜ¹¹ 
    Cᶜ¹² = l₀2 * Cᶜ¹² 
    Cᶜ¹³ = l₀2 * Cᶜ¹³ 
    Cᶜ¹⁴ = l₀2 * Cᶜ¹⁴ 
    Cᶜ²¹ = l₀2 * Cᶜ²¹ 
    Cᶜ²² = l₀2 * Cᶜ²² 
    Cᶜ²³ = l₀2 * Cᶜ²³ 
    Cᶜ²⁴ = l₀2 * Cᶜ²⁴ 
    Cᶜ³¹ = l₀2 * Cᶜ³¹ 
    Cᶜ³² = l₀2 * Cᶜ³² 
    Cᶜ³³ = l₀2 * Cᶜ³³ 
    Cᶜ³⁴ = l₀2 * Cᶜ³⁴ 
    Cᶜ⁴¹ = l₀2 * Cᶜ⁴¹ 
    Cᶜ⁴² = l₀2 * Cᶜ⁴² 
    Cᶜ⁴³ = l₀2 * Cᶜ⁴³ 
    Cᶜ⁴⁴ = l₀2 * Cᶜ⁴⁴ 
    
    Kᶜ¹¹ = l₀2 * Kᶜ¹¹
    Kᶜ¹² = l₀2 * Kᶜ¹²
    Kᶜ¹³ = l₀2 * Kᶜ¹³
    Kᶜ¹⁴ = l₀2 * Kᶜ¹⁴
    Kᶜ²¹ = l₀2 * Kᶜ²¹
    Kᶜ²² = l₀2 * Kᶜ²²
    Kᶜ²³ = l₀2 * Kᶜ²³
    Kᶜ²⁴ = l₀2 * Kᶜ²⁴
    Kᶜ³¹ = l₀2 * Kᶜ³¹
    Kᶜ³² = l₀2 * Kᶜ³²
    Kᶜ³³ = l₀2 * Kᶜ³³
    Kᶜ³⁴ = l₀2 * Kᶜ³⁴
    Kᶜ⁴¹ = l₀2 * Kᶜ⁴¹
    Kᶜ⁴² = l₀2 * Kᶜ⁴²
    Kᶜ⁴³ = l₀2 * Kᶜ⁴³
    Kᶜ⁴⁴ = l₀2 * Kᶜ⁴⁴
    
    contact_energy = l₀2 * contact_energy
    
    if exact
        
        Θ₁ = toangle(ΔR₁)
        Θ₂ = toangle(ΔR₂)
        
        Tₛ⁻¹Θ₁ = Tₛ⁻¹(Θ₁)
        Tₛ⁻¹Θ₂ = Tₛ⁻¹(Θ₂)
        
        Cᶜ²¹ = Cᶜ²¹ * Tₛ⁻¹Θ₁'
        Cᶜ²² = Cᶜ²² * Tₛ⁻¹Θ₁' 
        Cᶜ²³ = Cᶜ²³ * Tₛ⁻¹Θ₁'
        Cᶜ²⁴ = Cᶜ²⁴ * Tₛ⁻¹Θ₁'
        Cᶜ⁴¹ = Cᶜ⁴¹ * Tₛ⁻¹Θ₂'
        Cᶜ⁴² = Cᶜ⁴² * Tₛ⁻¹Θ₂' 
        Cᶜ⁴³ = Cᶜ⁴³ * Tₛ⁻¹Θ₂'
        Cᶜ⁴⁴ = Cᶜ⁴⁴ * Tₛ⁻¹Θ₂'
        
    end
    
    Tᶜ = [Tᶜ¹; Tᶜ²; Tᶜ³; Tᶜ⁴]
    Kᶜ = hcat(vcat(Kᶜ¹¹, Kᶜ²¹, Kᶜ³¹, Kᶜ⁴¹), vcat(Kᶜ¹², Kᶜ²², Kᶜ³², Kᶜ⁴²), vcat(Kᶜ¹³, Kᶜ²³, Kᶜ³³, Kᶜ⁴³), vcat(Kᶜ¹⁴, Kᶜ²⁴, Kᶜ³⁴, Kᶜ⁴⁴))
    Cᶜ = hcat(vcat(Cᶜ¹¹, Cᶜ²¹, Cᶜ³¹, Cᶜ⁴¹), vcat(Cᶜ¹², Cᶜ²², Cᶜ³², Cᶜ⁴²), vcat(Cᶜ¹³, Cᶜ²³, Cᶜ³³, Cᶜ⁴³), vcat(Cᶜ¹⁴, Cᶜ²⁴, Cᶜ³⁴, Cᶜ⁴⁴))
    
    return contact_energy, Tᶜ, Kᶜ,  Cᶜ
    
end

function assemble_contact!(conf::BeamsConfiguration, state::SimulationState, master_surface, slave_surface::BeamElementSurface, properties_contact, inter, params::SimulationParams) 
    
    # Unpack information 
    @unpack nᴳ_beams, ωᴳ_beams, zᴳ_beams, α = params
    gauss_params  = (nᴳ_beams, ωᴳ_beams, zᴳ_beams)   
    
    # Unpack configuration to retrieve node and beam information
    @unpack nodes, beams = conf
    
    # Unpack contact information
    @unpack contact_beams = slave_surface
    
    # Loop over all beam beams to compute their contributions
    for ind in contact_beams
        
        # Retrieve corresponding beams
        beam = beams[ind]
        
        # Retrieve node indices for the current beam beam
        n1 = beam.node1
        n2 = beam.node2
        
        # Extract nodal data for the two nodes connected by the beam
        X₁, X₂ = nodes.X₀[n1], nodes.X₀[n2]             # Initial positions
        u₁, u₂ = nodes.u[n1], nodes.u[n2]               # Displacements
        u̇₁, u̇₂ = nodes.u̇[n1], nodes.u̇[n2]           # Velocities
        ẇ₁, ẇ₂ = nodes.ẇ[n1], nodes.ẇ[n2]           # Rotational velocities
        R₁, R₂ = nodes.R[n1], nodes.R[n2]               # Rotational transformations
        ΔR₁, ΔR₂ = nodes.ΔR[n1], nodes.ΔR[n2]
        
        # Pack initialization constants for the beam computation
        init = (X₁, X₂, beam.l₀, beam.Rₑ⁰)  # Beam index, geometry, and initial rotation
        constants = (init, gauss_params, beam.properties.radius, properties_contact, master_surface)
        
        # Compute beam-level contributions
        contact_energy, Tᶜ, Kᶜ, Cᶜ = compute_contact_beams(u₁, u₂, R₁, R₂, ΔR₁, ΔR₂, u̇₁, u̇₂, ẇ₁, ẇ₂, constants) 
        
        # Assemble contributions into global vectors and matrices
        dofs1 = nodes.global_dofs[n1]               # dofs for node 1
        dofs2 = nodes.global_dofs[n2]               # dofs for node 2
        dofs = vcat(dofs1, dofs2)              # Combine dofs for the beam
        
        # Accumulate beam-level energies
        state.energyⁿ⁺¹.contact_energy += contact_energy
        
        # Add beam-level forces to global force vectors
        state.forcesⁿ⁺¹.Tᶜ[dofs] += Tᶜ               # Contact forces
        
        # Add beam-level matrices to global matrices
        state.matricesⁿ⁺¹.K.nzval[beam.global_sparsity_map] += vec(- Kᶜ)  # Stiffness matrix
        state.matricesⁿ⁺¹.C.nzval[beam.global_sparsity_map] += vec(- (1 + α) * Cᶜ)  # Damping matrix
        
    end
    
end 