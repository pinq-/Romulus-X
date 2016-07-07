include("leikkausvoimien_laskenta.jl")
include("leikkausvoimat_taivutus.jl")
include("lahto_arvot.jl")
include("apufunktiot.jl")
using Roots
using Gadfly
# Ratkaistaan leikkaustason kulma φ
#φ_1, φ_2 = func_φ(β_friction,α_rake,Z)

# Lasketaan leikkausvenymä γ
#γ = func_γ(α_rake,φ_1)

# Lasketaan horisontaalinen leikkausvoima
#Fc_1 = func_Fc(β_friction,φ_1,α_rake,τ,w,γ,t0,R)
#Fc_2 = func_Fc(β_friction,φ_2,α_rake,τ,w,γ,t0,R)


###############################
# φ:n laskennan tarkistusta
###############################
#φ:n ratkaiseminen numeerisesti (Ei tarvitse laskea, funktio func_φ toimii)
#f(φ) = (1.0 - (sin(β_friction)*sin(φ)) / (cos(β_friction-α_rake)*cos(φ-α_rake))) * (1.0/(cos(φ-α_rake)^2) - 1.0/(sin(φ)^2)) + ( cot(φ)+tan(φ-α_rake)+Z ) * ( (sin(β_friction)/(cos(β_friction-α_rake))) * ( cos(φ)/(cos(φ-α_rake)) + (sin(φ)*sin(φ-α_rake))/(cos(φ-α_rake)^2)) )
#fzeros(f,-2,2)

###############################
# φ:n piirto α_rake:a vastaan
###############################
#Lasketaan Z annetuilla materiaaliarvoilla ja lastun paksuudella
Z1 = func_Z(R,τ,t0)

# Määritellään funktiot, jotka piirretään
φ_1plot(α_rake) = acot( tan(β_friction - α_rake) + sqrt(1.0 + tan(β_friction - α_rake)^2 + Z1*(tan(β_friction - α_rake) + tan(α_rake))) )
φ_2plot(α_rake) = acot( tan(β_friction - α_rake) - sqrt(1.0 + tan(β_friction - α_rake)^2 + Z1*(tan(β_friction - α_rake) + tan(α_rake))) )

#Etsitään neliöjuuren sisällä olevan osan merkinvaihdot, jotta saadaan x-akselin rajat min1 & max1
n_juuri(α_rake) = 1.0 + tan(β_friction - α_rake)^2 + Z1*(tan(β_friction - α_rake) + tan(α_rake))
ii = 0
n_juuri_nollat = fzeros(n_juuri,-2,2)
while ii < 1
    n_juuri_nollat = fzeros(n_juuri,-2,2)
    if length(n_juuri_nollat) < 4
        n_juuri_nollat = fzeros(n_juuri,-2,2)
    elseif length(n_juuri_nollat) >= 4
        ii = 1
    end
end
#Etsitään Origon molemmilla puolilla olevat, lähimmät nollakohdat
j = 2
min1 = n_juuri_nollat[1]
max1 = n_juuri_nollat[length(n_juuri_nollat)]

while j <= length(n_juuri_nollat)
    min2 = n_juuri_nollat[j]
    if min2 < 0 && min2 >= min1
        min1 = min2
    elseif min2 >= 0 && min2 < max1
        max1 = min2
    end
    j = j +1
end

# Muuttuja jota käytetään varmistamaan, ettei pyöristysvirheet aiheuta ongelmia min1:ssä ja max1:ssä
jj = 0.00001
min1 = min1 + jj
max1 = max1 - jj

#Plottaus radiaaneissa
plot1 = plot(
    # Kaksi viivaa, joten erotellaan layereihin
    # first layer
    layer(
        # Määritellään funktio, sekä minimi- ja maksimirajat jolla välillä funktio piirretään
        φ_1plot,min1,max1,
        # Määritellään viivan väri
        Theme(
            default_color=colorant"orange")),
    # second layer
    layer(
        φ_2plot,min1,max1,
        Theme(
            default_color=colorant"purple")),
    # x- ja y-akselien otsikkotiedot, sekä määritellään x-akselin ala- ja yläraja
    Guide.ylabel("φ [rad]"),
    Guide.xlabel("α [rad]"),
    Scale.x_continuous(minvalue=-2, maxvalue=2)
)

################################
# Leikkausvoimien laskenta leikkaussyvyyden funktiona
################################
# Leikkaussyvyys
t1 = [0.00000005:0.00000005:0.000001]
# Alustetaan muuttujat
Fc = Array(Float64,t1)
Fc_w = Array(Float64,Fc)
Fv = Array(Float64,Fc)
Fs = Array(Float64,Fc)
i = 1

while i  <= length(t1)
    # Päivitetään leikkaussyvyys
    local t0::Float64
    t0 = t1[i]
    # Lasketaan Z
    local Z::Float64
    Z = func_Z(R,τ,t0)
    # Ratkaistaan leikkaustason kulma φ
    φ_1, φ_2 = func_φ(β_friction,α_rake,Z)
    # Lasketaan leikkausvenymä
    local γ::Float64
    γ = func_γ(α_rake,φ_1)
    # Lasketaan leikkausvoimat
    Fc[i] = func_Fc(β_friction,φ_1,α_rake,τ,w,γ,t0,R)
    Fc_w[i] = func_Fc_w(Fc[i],w)
    Fv[i] = func_Fv(Fc[i],β_friction,α_rake)
    Fs[i] = func_Fs(Fc[i],φ_1,Fv[i])
    # Päivitetään laskuri ennen seuraavaa kierrosta
    i = i + 1
end

################################
# Leikkausvoiman Fc plottaus leikkaussyvyyden funktiona
################################
plot2 = plot(
    # x:n arvot leikkaussyvyydeksi, ja muutetaan [mm]
    x = t1.*1000,
    # y:n arvot leikkausvoimaksi per terä mm [N/mm]
    y=Fc_w,
    # Määritellään, että piirretään viiva, eikä datapisteitä
    Geom.line,
    # x- ja y-akselien otsikkotiedot
    Guide.xlabel("t0 [mm]"),
    Guide.ylabel("Fc/w [N/mm]"),)


###############################
# Kitkakulman plottaus leikkausnopeuden funktiona
###############################
plot3 = plot(
    β_friction_func,1,1000,
    Scale.x_log10,
    Scale.x_continuous(minvalue=1, maxvalue=1000),
    # x- ja y-akselien otsikkotiedot
    Guide.xlabel("v [mm/s]"),
    Guide.ylabel("β [deg]"),
    )


###############################
# Leikkausvoiman laskenta ottaen huomioon terän tylsyyden
###############################

Z = func_Z(R,τ,t0)
# Lasketaan stagnant materiaalialueen pituus terän edessä
L_AB = func_L_AB(r,α_rake,λ) #[m]
# Lasketaan slip-line kulma
ξ = func_ξ(m_plough) #[rad]
# Lasketaan prow-kulma
ρ = func_ρ(λ,ξ) #[rad]
# Lasketaan työkalun kärjen ja työkappaleen välinen kitkakulma
β_friction_plough = func_β_fric_plough(ξ,ρ,λ) #[rad]
#=
φ_plough = fzeros(func_φ_plough,0,6)

if length(φ_plough) == 0
    println("φ_plough ei löytynyt nollakohtia!")
elseif length(φ_plough) > 1
        φ_plough = φ_plough[1]
end
=#
Fc_plough = Array(Float64,t1)
Fc_plough_w = Array(Float64,Fc)
i2 = 1
while i2  <= length(t1)
    # Päivitetään leikkaussyvyys
    global t0_global
    t0_global = t1[i2]
    # Lasketaan Z
    global Z_global
    Z_global = func_Z(R,τ,t0_global)
    # Ratkaistaan leikkaustason kulma φ
    φ_plough = fzeros(func_φ_plough,0,4)
    if length(φ_plough) == 0
        println("φ_plough ei löytynyt nollakohtia!")
    elseif length(φ_plough) > 1
            φ_plough = φ_plough[1]
    end
    # Lasketaan leikkausvenymä
    γ = func_γ(α_rake,φ_plough[1])
    # Lasketaan Q
    Q = func_Q(β_friction,φ_plough[1],α_rake)
    # Lasketaan leikkausvoimat
    Fc_plough[i2] = func_Fc_plough(τ,w,γ,α_rake,t0_global,φ_plough[1],Q,m_plough,L_AB,ξ,β_friction_plough,λ,β_friction,R)
    Fc_plough_w[i2] = func_Fc_w(Fc_plough[i2],w)
    # Päivitetään laskuri ennen seuraavaa kierrosta
    i2 = i2 + 1
end

################################
# Leikkausvoiman Fc plottaus leikkaussyvyyden funktiona
################################
plot4 = plot(
    # first layer
    layer(
        # x:n arvot leikkaussyvyydeksi, ja muutetaan [mm]
        x = t1.*1000,
        # y:n arvot leikkausvoimaksi per terä mm [N/mm]
        y=Fc_plough_w,
        # Määritellään viivan väri
        Theme(
            default_color=colorant"blue"),
        # Määritellään, että piirretään viiva, eikä datapisteitä
        Geom.line,
        ),
    # second layer
    layer(
        # x:n arvot leikkaussyvyydeksi, ja muutetaan [mm]
        x = t1.*1000,
        # y:n arvot leikkausvoimaksi per terä mm [N/mm]
        y=Fc_w,
        Theme(
            default_color=colorant"red"),
        # Määritellään, että piirretään viiva, eikä datapisteitä
        Geom.line,
        ),
    # Määritellään, että piirretään viiva, eikä datapisteitä
    #Geom.line,
    # x- ja y-akselien otsikkotiedot
    Guide.xlabel("t0 [mm]"),
    Guide.ylabel("Fc_plough/w [N/mm]"),)


################################
# Leikkausvoiman Fc laskenta ja plottaus kun sorvataan ideaalia pölliä 2D-tapaus
################################

sorvaus(R_alku,R_purilas,t0)

#=
#Terän etäisyys globaalista origosta (Karan keskipiste)
L_terä = [0.5:-t0:0.0] #[m]

Θ_alku = (R_alku/t0) * 2*pi
Θ_temp = [Θ_alku:-0.01:0.0]
i3 = 1

t0_temp = zeros(length(Θ_temp))
Θ_kuljettu = zeros(length(Θ_temp))
R_temp = zeros(length(Θ_temp))
Fc_sorvaus = zeros(length(Θ_temp))
Fc_w_sorvaus = zeros(length(Θ_temp))
Fv_sorvaus = zeros(length(Θ_temp))
Fs_sorvaus = zeros(length(Θ_temp))

for L_terä in [0.5:-t0:0.0]
    if L_terä >= R_pölli
        Fc[] = 0.0
    elseif L_terä <= R_pölli
        # Lasketaan leikkaussyvyys t0
        t0_temp = R_pölli - L_terä

        # Lasketaan Z
        local Z
        Z = func_Z(R,τ,t0_temp[i3])
        # Ratkaistaan leikkaustason kulma φ
        φ_1, φ_2 = func_φ(β_friction,α_rake,Z)
        # Lasketaan leikkausvenymä
        local γ
        γ = func_γ(α_rake,φ_1)
        # Lasketaan leikkausvoimat
        Fc_sorvaus[i3] = func_Fc(β_friction,φ_1,α_rake,τ,w,γ,t0_temp[i3],R)
        Fc_w_sorvaus[i3] = func_Fc_w(Fc_sorvaus[i3],w)
        Fv_sorvaus[i3] = func_Fv(Fc_sorvaus[i3],β_friction,α_rake)
        Fs_sorvaus[i3] = func_Fs(Fc_sorvaus[i3],φ_1,Fv_sorvaus[i3])

        i3 = i3 + 1
=#
#=
# Looppi laskee leikkausvoimat kuljetun kulman Θ_kuljettu avulla
while i3 <= length(Θ_temp)
    Θ_kuljettu[i3] = Θ_alku - Θ_temp[i3]
    R_temp[i3] = R_alku - t0 / (2*pi) * Θ_kuljettu[i3]

    if R_temp[i3] <= R_purilas
        break
    else
        # Päivitetään leikkaussyvyys
        if Θ_kuljettu[i3] <= 2*pi
            t0_temp[i3] = R_alku - R_temp[i3]
            if t0_temp[i3] <= 0
                # Lasketaan leikkausvoimat
                Fc_sorvaus[i3] = 0.0
                Fc_w_sorvaus[i3] = 0.0
                Fv_sorvaus[i3] = 0.0
                Fs_sorvaus[i3] = 0.0

                i3 = i3 + 1
                continue
            end
        else
            t0_temp[i3] = t0
        end
        # Lasketaan Z
        local Z
        Z = func_Z(R,τ,t0_temp[i3])
        # Ratkaistaan leikkaustason kulma φ
        φ_1, φ_2 = func_φ(β_friction,α_rake,Z)
        # Lasketaan leikkausvenymä
        local γ
        γ = func_γ(α_rake,φ_1)
        # Lasketaan leikkausvoimat
        Fc_sorvaus[i3] = func_Fc(β_friction,φ_1,α_rake,τ,w,γ,t0_temp[i3],R)
        Fc_w_sorvaus[i3] = func_Fc_w(Fc_sorvaus[i3],w)
        Fv_sorvaus[i3] = func_Fv(Fc_sorvaus[i3],β_friction,α_rake)
        Fs_sorvaus[i3] = func_Fs(Fc_sorvaus[i3],φ_1,Fv_sorvaus[i3])

        i3 = i3 + 1
    end
end
=#

plot5 = plot(
    # first layer
    layer(
        # x:n arvot leikkaussyvyydeksi, ja muutetaan [mm]
        x = t0_temp[1:i3].*1000,
        # y:n arvot leikkausvoimaksi per terä mm [N/mm]
        y=Fc_w_sorvaus[1:i3],
        # Määritellään viivan väri
        Theme(
            default_color=colorant"blue"),
        # Määritellään, että piirretään viiva, eikä datapisteitä
        Geom.line,
        ),
    # x- ja y-akselien otsikkotiedot
    Guide.xlabel("t0 [mm]"),
    Guide.ylabel("Fc_sorvaus/w [N/mm]"),)

plot6 = plot(
    # first layer
    layer(
        # x:n arvot [rad]
        x = Θ_kuljettu[1:i3],
        # y:n arvot leikkausvoimaksi per terä mm [N/mm]
        y = Fc_w_sorvaus[1:i3],
        # Määritellään viivan väri
        Theme(
            default_color=colorant"blue"),
        # Määritellään, että piirretään viiva, eikä datapisteitä
        Geom.line,
        ),
    # x- ja y-akselien otsikkotiedot
    Guide.xlabel("Θ_kuljettu [rad]"),
    Guide.ylabel("Fc_sorvaus/w [N/mm]"),)

#=
################################
# Leikkausvoiman Fc laskenta taivutustapauksessa
################################
# Lasketaan materiaalin plane strain modulus
E_pilkku = func_E_pilkku(E,ν)
# Lasketaan käyristymissäde jolla leikkaus muuttuu elastis-plastiseksi
Rp = func_Rp(t0,E,σy)

# Lasketaan kiertymä murtuman juuressa
Θ_root = func_Θ_root(E,t0,χ,R,δ_cut)
# Lasketaan käyristymisssäde murtuman juuressa
R_root = func_R_root(χ,t0,Θ_root)

# Tarkistetaan onko elastinen vai elastis-plastinen Leikkaus
# Jos R_root <= Rp => elastis-plastinen == TRUE
tapaus = func_el_pl_case(R_root,Rp)


if tapaus == 0 #Elastinen leikkaus
    println("Elastinen leikkaus")
    X = func_X(E,t0,χ,R)
    GtGc = func_GtGc(X,δ_cut)
    # Tarkistetaan tapahtuuko leikkauksessa teräkosketusta
    if func_teräkosketus_elast(GtGc,δ_cut) == 0 #Ei teräkosketusta
        println("Ei teräkosketusta")
        # Lasketaan leikkausvoima Fc
        Fc = func_Fc_bend_e(R,w,Ga,α_rake,β_friction)
    elseif func_teräkosketus_elast(GtGc,δ_cut) == 1 #Teräkosketus
        println("Teräkosketus")
        # Lasketaan leikkausvoima Fc
        Fc = func_Fc_crack_e(R,w,Ga,β_friction,E_pilkku,t0,χ,α_rake)
    end
elseif tapaus == 1 #Elastis-plastinen leikkaus
    println("Elastis-plastinen leikkaus")
    # Lasketaan säteiden suhde k0
    k0 = func_k0(Rp,R_root)
    # Tarkistetaan tapahtuuko leikkauksessa teräkosketusta
    if func_teräkosketus_plast(χ,σy,δ_cut,E,k0) == 0 #Ei teräkosketusta
        println("Ei teräkosketusta")

        # Lasketaan apukerroin γ1_williams
        γ1_williams = func_γ1_williams(k0)
        # Lasketaan G_hattu
        G_hattu = func_G_hattu(σy,t0,E)

        # Lasketaan leikkausvoima Fc
        Fc = func_Fc_plastic(Ga,α_rake,β_friction,G_hattu,γ1_williams,R,w)
    elseif func_teräkosketus_plast(χ,σy,δ_cut,E,k0) == 1 #Teräkosketus
        println("Teräkosketus")

        # Lasketaan G_hattu
        G_hattu = func_G_hattu(σy,t0,E)

        # Lasketaan leikkausvoima Fc
        Fc = func_Fc_crack_p(R,w,G_hattu,k0)
    end
end
=#








######################
#EOF
