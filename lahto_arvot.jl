include("apufunktiot.jl")
# SI-yksiköt
################################
#          Terägeometria       #
################################
#teroituskulma
const β_edge_deg = 20.0                #[deg]
const β_edge = func_radians(β_edge_deg)               #[rad]

#Päästökulma (clearance angle) (Muuttuu simulaation aikana, tarvitaan tieto sorvimallista)
α_clearance_deg = 1.0                           #[deg]
α_clearance = func_radians(α_clearance_deg)     #[rad]

#Leikkauskulma == Terän rintapinnan kulma leikkauspinnasta
δ_cut = β_edge + α_clearance                #[rad]

# Rake angle
α_rake = pi/2 - δ_cut                       #[rad]

# Leikkauspituus (Terän pituus, joka koskettaa pölliä)
const w = 1.350                                     #[m]

# Terän kärjen pyöristyssäde
const r = 5.0e-6                                  #[m]

################################
#        Materiaaliarvot       #
################################
# Shear stress along the shear plane
const τ = 470.0e6                      #[N/m^2]
# Fracture toughness
const R = 5.0e3                        #[J/m^2] = [Nm]

# Kimmomoduli (Tangentiaalisessa suunnassa)
const E = 0.113 * 9.7e9               #[N/m^2]

# Poissonin vakio
const ν = 0.362                       #[-]

# Yield stress
const σy = 80e3                       #[N/m^2]

# Pöllin säde
const R_alku = 0.120                  #[m]

# Purilaan säde
const R_purilas = 0.105               #[m]

#=
# Kitkakerroin puun ja terän välillä (Muuttuu leikkausnopeuden funktiona semi-log skaalassa, MUTTA leikkausvoima ei vaikuta merkittävästi!!!!)
μ = 0.3                         #[-]

# Angle of friction
β_friction = atan(μ)            #[rad]
=#
################################
# Angle of frictionin muutos leikkausnopeuden funktiona
const b = 28.0
const m = -4.6666666666667
β_friction_func(v) = m*log10(v) + b #[deg]

################################
#        Leikkausarvot         #
################################
# Leikkaamattoman lastun paksuus (Viilun paksuus jonka sorvin käyttäjä valitsee)
t0 = 1.5                          #[mm]
#muutetaan SI-yksiköihin
t0 = t0/1000                        #[m]

# Leikkausnopeus (Nopeus, jolla terä liikkuu materiaalia kohti leikkaussuunnassa,
# muuttuu simulaation aikana, pitää laskea pöllin kehänopeudesta ja terän lineaariliikkeestä)
v = 180                            #[m/s]

# Muutetaan mm/s kitkakulman laskemista varten
v_mm = v * 1000                 #[mm/s]
if v_mm == 0
    β_friction = 0
else
    β_friction = func_radians(β_friction_func(v_mm)) #[rad]
end

#Ylikirjotetaan arvoja testejä varten
#β_friction = func_radians(37)
#α_rake = func_radians(35)

#Terän pyöristyssäteen laskentaan
const λ = 0.0
const m_plough = 0.099 #Kitkakerroin, arvot 0...1.0, paitsi että φ_plough laskenta ei toimi jos arvo on 1

#Taipuisan palkin laskentaan
const χ = 0.64
const Ga = 0














######################
#EOF
