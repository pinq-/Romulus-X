include("leikkausvoimien_laskenta.jl")
include("lahto_arvot.jl")
include("apufunktiot.jl")

using PyPlot

#Animaatiofunktio
function sorvaus(R_alku::Float64,R_purilas::Float64,t0::Float64)
    # Kuvataan profiili polaarikoordinaateissa. Jokainen piste on muotoa (Θ,r)
    # Jaetaan profiili n+1:een pisteeseen
    # Pisteiden lukumäärä
    n::Int64 = 359
    # Yhden sektorin kulma
    #Θ_sektori = 2*pi / n
    R_sektori = Array(Float64, n) #zeros(n)
    Θ_sektori = Array(Float64, n) #zeros(n)
    # Jaetaan profiili n:ään pisteeseen
    @simd for i in range(1,n)
        @fastmath @inbounds R_sektori[i] = R_alku * (rand(50:115)/100)
        @fastmath @inbounds Θ_sektori[i] = (2*pi / n) * i
    end
    # Lisätään profiilin dataan alkupiste, jotta plottauksessa profiili ei ole aukinainen
    unshift!(Θ_sektori,0.0)
    unshift!(R_sektori,R_sektori[n])


    fig = plt[:figure]()
    ax = fig[:add_subplot](111, projection="polar")
    # Määritellään polaarikoordinaatiston radiaalinen raja
    ax[:set_ylim](0,R_alku+0.03)

    # Kulma, jonka verran pöllin profiili pyörii per askel
    ω::Float64 = func_radians(1.1) #[rad]

    # Terän paikka koordinaatistossa
    tera_asema::Float64 = 0.123 #R_sektori[n]
    tera_kulma::Float64 = 0.0 #[rad]

    # Plotataan profiili, lasketaan joka aika-askeleella pisteille uusi kulma-asema ja päivitetään plottaus
    i::Int64 = 1
    Fc_sorvaus = Float64[]

    while R_sektori[1] >= R_purilas #tera_asema >= R_purilas
        if i == 1
            fig[:canvas][:draw]()
            # Piirretään pöllin alkuasema
            profile, = ax[:plot](Θ_sektori,R_sektori)
            # Terän alkuasema
            tera, = ax[:plot](tera_kulma,tera_asema,"k.")
            # Alkupiste pölliin, niin näkee pyörimisen paremmin
            #profiili_alku, = ax[:plot](Θ_sektori[1],R_sektori[1],"k.")

        elseif i == 2
            #Fc_sorvaus_i = 0
            #push!(Fc_sorvaus,Fc_sorvaus_i)

            # Lasketaan pöllin uusi kulma-asema
            #global Θ_sektori_uus
            @fastmath Θ_sektori_uus = Θ_sektori.-(ω*i)
            # Lasketaan terän uusi asema
            @fastmath tera_asema_uus::Float64 = tera_asema - (t0 / (2*pi) * ω*i)
            # Asetetaan datapisteille uudet asemat
            profile[:set_data](Θ_sektori_uus,R_sektori)
            #profiili_alku[:set_data](Θ_sektori_uus[1],R_sektori[1])
            tera[:set_data](tera_kulma,tera_asema_uus)

        else
            # Katotaan kuinka moni piste menee terälinjan ohi per askel
            piste = Int64[]
            @fastmath Θ_sektori_edellinen = Θ_sektori_uus.-ω
            @simd for k in range(1,n+1)
                #global a1,a2
                @fastmath @inbounds a1 = func_cartesis(R_sektori[k],Θ_sektori_edellinen[k]-tera_kulma)
                @fastmath @inbounds a2 = func_cartesis(R_sektori[k],Θ_sektori_uus[k]-tera_kulma)
                if sign(a1[1]) == sign(a2[1]) == 1 && sign(a1[2]) != sign(a2[2])
                    push!(piste,k)
                end
            end
            #Tsekataan kulkiko yksikään piste terän ohi
            #if length(piste) == 0

            # Ehto joka tsekkaa onko terä pöllin sisäpuolella IF TRUE -> lasketaan leikkausvoima ja lasketaan pöllin pisteelle uusi koordinaatti
            @simd for l in range(1,length(piste))
                if tera_asema_uus <= R_sektori[piste[l]]
                    # Lasketaan leikkauspaksuus
                    local t0_temp::Float64
                    @fastmath t0_temp = R_sektori[piste[l]] - tera_asema_uus
                    # Lasketaan Z
                    local Z::Float64
                    @fastmath Z = func_Z(R,τ,t0_temp)
                    # Ratkaistaan leikkaustason kulma φ
                    local φ_1::Float64, φ_2::Float64
                    @fastmath φ_1, φ_2 = func_φ(β_friction,α_rake,Z)
                    # Lasketaan leikkausvenymä
                    local γ::Float64
                    @fastmath γ = func_γ(α_rake,φ_1)
                    # Lasketaan leikkausvoimat
                    local Fc_sorvaus_i::Float64
                    @fastmath Fc_sorvaus_i = func_Fc(β_friction,φ_1,α_rake,τ,w,γ,t0_temp,R)
                    push!(Fc_sorvaus,Fc_sorvaus_i)
                    #Fc_w_sorvaus[i] = func_Fc_w(Fc_sorvaus[i],w)
                    #Fv_sorvaus[i] = func_Fv(Fc_sorvaus[i],β_friction,α_rake)
                    #Fs_sorvaus[i] = func_Fs(Fc_sorvaus[i],φ_1,Fv_sorvaus[i])

                    # Lasketaan leikatun pisteen uusi asema
                    @fastmath R_sektori[piste[l]] = R_sektori[piste[l]] - t0_temp
                else
                    Fc_sorvaus_i = 0
                    push!(Fc_sorvaus,Fc_sorvaus_i)
                end
                # Lasketaan pöllin uusi kulma-asema
                #global Θ_sektori_uus
                @fastmath Θ_sektori_uus = Θ_sektori.-(ω*i)
                # Lasketaan terän uusi asema
                @fastmath tera_asema_uus = tera_asema - (t0 / (2*pi) * ω*i)
                # Asetetaan datapisteille uudet asemat
                profile[:set_data](Θ_sektori_uus,R_sektori)
                #profiili_alku[:set_data](Θ_sektori_uus[1],R_sektori[1])
                tera[:set_data](tera_kulma,tera_asema_uus)
            end
        end
        i = i + 1

        # Pausetetaan animaatiota, jotta sen näkee
        #plt[:pause](0.01) #[s] #Matplotlib:n sleep komento, hiton hidas, minimi pause on 0.01 s.
        # Matplotlib:n pause komento eriteltynä animaation nopeuttamiseksi, piirretään uusiksi vain tarvittavat osat.
        #ax[:draw_artist](profile)
        #ax[:draw_artist](profiili_alku)
        #ax[:draw_artist](tera)
        fig[:canvas][:update]()
        fig[:canvas][:flush_events]()
        #sleep(0.002) #Julian oma sleep komento. Minimiaika on 1 ms

    end
end #Funktion end

sorvaus(R_alku,R_purilas,t0)
#EOF
