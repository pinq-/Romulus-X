# Atkinsin muokkaama Ernst-Merchantin leikkausteoria. Lastu leikkautuu irti kappaleesta leikkaustason kohdalta.
# Täysin terävä terä.

#Laskee termin Z [-]
function func_Z(R::Float64,τ::Float64,t0::Float64)
    Z::Float64 = R / (τ * t0)
end

# Laskee leikkaustason kulman φ [rad]
function func_φ(β_friction::Float64,α_rake::Float64,Z::Float64)
    φ_1::Float64 = acot( tan(β_friction - α_rake) + sqrt(1.0 + tan(β_friction - α_rake)^2 + Z*(tan(β_friction - α_rake) + tan(α_rake))) )
    φ_2::Float64 = acot( tan(β_friction - α_rake) - sqrt(1.0 + tan(β_friction - α_rake)^2 + Z*(tan(β_friction - α_rake) + tan(α_rake))) )

    #println("φ_1 = ",φ_1*180/pi, " [deg]")
    #println("φ_2 = ",φ_2*180/pi, " [deg]")

    return φ_1, φ_2
    #φ = jotain
end

#Laskee leikkausvenymän γ leikkaustasossa φ
function func_γ(α_rake::Float64,φ::Float64)
    γ::Float64 = cos(α_rake) / ( sin(φ)*cos(φ - α_rake) )
end

# Laskee leikkauspinnan suuntaisen leikkausvoiman komponentin Fc
function func_Fc(β_friction::Float64,φ::Float64,α_rake::Float64,τ::Float64,w::Float64,γ::Float64,t0::Float64,R::Float64)
    # Q = friction correction
    Q::Float64 = 1.0 - ( sin(β_friction)*sin(φ) ) / ( cos(β_friction - α_rake) * cos(φ - α_rake) )
    Fc::Float64 = (τ*w*γ / Q) * t0 + (R*w / Q) #[N]
    #println("Fc = ",Fc," [N]")
    return Fc
end

# Laskee leikkauspintaan kohtisuoraan suuntaisen leikkausvoiman komponentin Fv
function func_Fv(Fc::Float64,β_friction::Float64,α_rake::Float64)
    Fv::Float64 = Fc * tan(β_friction - α_rake) #[N]
end

# Laskee horisontaalisen leikkausvoiman per leikkauspituus ja muuntaa SI-yksiköistä mm:hin
function func_Fc_w(Fc::Float64,w::Float64)
    Fc_w::Float64 = Fc / (1000.0 * w) #[N/mm]
end

# Laskee leikkausvoiman Fs joka vaikuttaa leikkaustasossa φ
function func_Fs(Fc::Float64,φ::Float64,Fv::Float64)
    Fs::Float64 = Fc*cos(φ) - Fv*sin(φ) #[N]
end


##############################################
# Atkinsin malli, johon on lisätty terän pyöristyssäteen vaikutus
##############################################

# Laskee horisontaalisen leikkausvoiman, jossa otetaan huomioon terän pyöristys
function func_Fc_plough(τ,w,γ,α_rake,t0,φ,Q,m_plough,L_AB,ξ,β_friction_plough,λ,β_friction,R) #[N]
    A = τ*w*γ*t0 / Q
    B = m_plough*τ*L_AB*w / Q
    C = sin(ξ - λ) / sin(ξ)
    D = sin(β_friction_plough + λ)*sin(β_friction)*sin(φ) / (sin(β_friction_plough)*cos(β_friction - α_rake)*cos(φ - α_rake))
    E = R*w/Q
    Fc_plough = A + B * ( C - D ) + E
end

# Laskee leikkauksen stagnant zonen pituuden
function func_L_AB(r,α_rake,λ)
    L_AB = r*(1+sin(α_rake)) / (cos(λ-α_rake)) #[m]
end

# Laskee terän ploughing voiman Fp kitkakulman β_friction_plough [rad]
function func_β_fric_plough(ξ,ρ,λ)
    β_friction_plough = atan( cos(2*ξ) / (1 + pi/2 - 2*ρ + 2*ξ - 2*λ + sin(2*ξ)) )
end

# Laskee prow-kulman ρ
function func_ρ(λ,ξ)
    ρ = asin( sin(λ) / (sqrt(2)*sin(ξ)) ) #[rad]
end

# Laskee slip-line kulman ξ
function func_ξ(m_plough)
    ξ = (1/2) * acos(m_plough) #[rad]
end

#Laskee kitkakertoimen Q
function func_Q(β_friction,φ,α_rake)
    Q = 1.0 - ( sin(β_friction)*sin(φ) ) / ( cos(β_friction - α_rake) * cos(φ - α_rake) )
end

# Leikkaustason kulma φ:tä ei voida ratkaista samalla kaavalla kuin edellä.
# Nyt tarvitsee ratkaista numeerisesti seuraavasta yhtälöstä
function func_φ_plough(φ_plough)
    A = sin(β_friction) / cos(β_friction - α_rake)
    B = cos(α_rake) / (sin(φ_plough)*cos(φ_plough - α_rake))
    C = m_plough*r*(1+sin(α_rake)) / (t0_global*cos(λ-α_rake))
    D = sin(ξ-λ) / sin(ξ)
    E = sin(β_friction_plough + λ)*sin(β_friction)*sin(φ_plough) / (sin(β_friction_plough)*cos(β_friction-α_rake)*cos(φ_plough-α_rake))

    F = sin(β_friction)*sin(φ_plough) / ( cos(β_friction - α_rake)*cos(φ_plough - α_rake) )
    G = cos(2*φ_plough - α_rake) / (sin(φ_plough)^2)
    H = m_plough*r*(1 + sin(α_rake)) / ( t0_global*cos(λ-α_rake) )
    I = sin(β_friction_plough + λ) * sin(β_friction) / ( sin(β_friction_plough)*cos(β_friction - α_rake) )

    φ_plough = A * ( B + C * (D - E) + Z_global ) - (1 - F)*( G + H*I ) #[rad]

    #φ_plough = sin(β_friction)/cos(β_friction - α_rake) * ( cos(α_rake)/(sin(φ_plough)*cos(φ_plough - α_rake)) + m_plough*r*(1+sin(α_rake))/(t0*cos(λ-α_rake)) * (sin(ξ-λ)/sin(ξ) - sin(β_friction_plough + λ)*sin(β_friction)*sin(φ_plough)/(sin(β_friction_plough)*cos(β_friction-α_rake)*cos(φ_plough-α_rake))) + Z ) - (1 - sin(β_friction)*sin(φ_plough) / ( cos(β_friction - α_rake)*cos(φ_plough - α_rake) ))*( cos(2*φ_plough - α_rake) / (sin(φ_plough)^2) + m_plough*r*(1 + sin(α_rake)) / ( t0*cos(λ-α_rake) )*sin(β_friction_plough + λ) * sin(β_friction) / ( sin(β_friction_plough)*cos(β_friction - α_rake) ) ) #[rad]
end
# Käytä tätä Roots-paketin funktiota laskemaan func_φ_plough:n nollakohdat
#fzeros(func_φ_plough,-2,2)

# Laskee lastun efektiivisen paksuuden, silloin kun λ != 0
function func_t_effective(t0,L_AB,λ)
    t_effective = t0 - L_AB*sin(λ) #[m]
end

# Laskee efektiivisen rake-kulman tapauksessa, jossa lastu on kontaktissa
# terän pyöristyksen kanssa, eikä rintapinnan kanssa. Riippuu Leikkaamattoman
# lastun paksuudesta
# EI OIKEIN TÄLLÄ HETKELLÄ!!! Jos ei ratkee, niin oleta λ == 0 jolloin tätä ja t_effective ei tarvita
function func_α_rake_effective(α_rake,r,t_effective)
    α_rake_effective = α_rake * sin( (r - t_effective) / r ) #[rad]
end










######################
#EOF
