include("leikkausvoimien_laskenta.jl")
include("3D_laskentafunktiot.jl")
include("lahto_arvot.jl")
include("apufunktiot.jl")
include("pollin_dynamics_funktiot.jl")
include("pollidata.jl")

using PyCall
using PyPlot
@pyimport matplotlib.cm as cm

#Animaatiofunktio
function sorvaus(R_alku::Float64,R_purilas::Float64,t0::Float64, R_sektori)
    # Kulma, jonka verran pöllin profiili pyörii per askel
    ω::Float64 = func_radians(11.1) #[rad]

    # Jaetaan pölli k:aan profiiliin
    k::Int64 = 20
    Z_pölli = func_z_serial(k,w)

    # Kuvataan profiili polaarikoordinaateissa. Jokainen piste on muotoa (Θ,r)
    # Jaetaan profiili n+1:een pisteeseen
    # Pisteiden lukumäärä
    n::Int64 = 359
    Θ_sektori, R_sektori_turha = func_sektorit_serial(n,k,R_alku)

    # Alustetaan matriisit yms.
    Θ_sektori_uus = copy(Θ_sektori)
    Fc_sorvaus_kootut = Array(Array{Float64},k-1)
    Fv_sorvaus_kootut = Array(Array{Float64},k-1)
    tempvector = Array(Float64,0)
    for i in range(1,k-1)
        Fc_sorvaus_kootut[i] = copy(tempvector)
        Fv_sorvaus_kootut[i] = copy(tempvector)
    end

    # Muutetaan profiilien pisteet karteesiseen koordinaatistoon
    X_pölli = Array(Float64, n+1, k)
    Y_pölli = Array(Float64, n+1, k)
    func_XY_pölli_serial!(X_pölli,Y_pölli,R_sektori, Θ_sektori, n, k)
    #println("X_pölli[1,1] = ", X_pölli[1,1])
    #println("Y_pölli[1,1] = ", Y_pölli[1,1])

    ##################################
    # Pöllin dynamiikan alustus
    ##################################
    # Laskee kaikkien profiilien pinta-alat
    A = func_area_elements(R_sektori, Θ_sektori, n, k)
    # Laskee kaikkien elementtien tilavuudet ja tallentaa ne matriisiin
    volume_elements = func_volume_elements(Z_pölli, A, n, k)
    # Laskee koko pöllin tilavuuden
    volume_total = func_volume_total(volume_elements)
    # laskee kaikkien elementtien tilavuuskeskiöt ja tallentaa ne Array of Arrays:iin
    CoV_elements = func_CoV_elements(X_pölli, Y_pölli, Z_pölli, A, n, k)
    # Laskee koko pöllin tilavuuskeskiön
    CoV_pölli = func_CoV_pölli(CoV_elements, volume_elements, volume_total, n, k)
    # Laskee kaikkien elementtien inertiat
    I_elements = func_inertias_elements(CoV_elements, X_pölli, Y_pölli, Z_pölli, ρ, n, k)
    # Laskee koko pöllin inertian
    I_total = func_inertias_pölli(CoV_elements, I_elements, volume_elements, ρ, CoV_pölli[1], CoV_pölli[2], CoV_pölli[3], n, k)
    ##################################

    # Terän paikka koordinaatistossa
    tera_asema = Array(Float64, k)
    Z_terä = Array(Float64, length(tera_asema))
    w_terä = w / (length(tera_asema)-1)
    #Terän päätepisteet. Jos eri niin looppi laskee terän pisteille oikean etäisyyden origosta
    tera_asema[1] = R_sektori[1,1]+0.01 #0.13
    tera_asema[k] = R_sektori[1,1]+0.01 #0.13
    tera_apukulma = atan( (tera_asema[k] - tera_asema[1])/w )
    # Lasketaan terän pituus
    @simd for i in range(1,length(tera_asema))
        if i == 1
            Z_terä[i] = 0.0
        elseif i < length(tera_asema)
            @fastmath @inbounds Z_terä[i] = w_terä*(i-1)
            @fastmath @inbounds tera_asema[i] = tera_asema[1] + (Z_terä[i] * tan(tera_apukulma))
        else
            @fastmath @inbounds Z_terä[i] = w_terä*(i-1)
        end
    end
    tera_kulma::Float64 = func_radians(0.0) #[rad]

    # Muutetaan terän pisteet karteesiseen koordinaatistoon
    X_terä = Array(Float64, length(tera_asema))
    Y_terä = Array(Float64, length(tera_asema))
    @simd for j in range(1,length(tera_asema))
            @fastmath @inbounds X_terä[j,1],Y_terä[j,1] = func_cartesis(tera_asema[j,1],tera_kulma)
    end



    # Plotataan profiili, lasketaan joka aika-askeleella pisteille uusi kulma-asema ja päivitetään plottaus
    ii::Int64 = 1

    # Määrittää kuvaajan
    fig = plt[:figure]()
    ax = fig[:add_subplot](111, projection="3d")
    # Määritellään koordinaatiston rajat
    ax[:set_xlim3d](-R_alku-0.53,R_alku+0.53)
    ax[:set_zlim3d](0,w+0.1)
    ax[:set_ylim3d](-R_alku-0.53,R_alku+0.53)

    while R_sektori[1,1] >= R_purilas
        if ii == 1
            ax[:view_init](elev=90, azim=-90)
            fig[:canvas][:draw]()
            # Piirretään pöllin alkuasema
            if k == 1
                profile = ax[:plot_wireframe](X_pölli, Y_pölli, Z_pölli)
            else
                profile = ax[:plot_surface](X_pölli, Y_pölli, Z_pölli)
            end
            # Terän alkuasema
            tera = ax[:plot_wireframe](X_terä,Y_terä, Z_terä, color="r")

        elseif ii == 2
            # Lasketaan pöllin uusi kulma-asema
            func_Θ_sektori_uus_serial!(Θ_sektori_uus,Θ_sektori,ω,ii)
            #println("ii = ", ii)
            #println("Θ_sektori_uus[1,1] = ", Θ_sektori_uus[1,1])

            # Muutetaan uus asema karteesiseen koordinaatistoon
            @simd for j in range(1,k)
                @simd for i in 1:(n+1)
                    @fastmath @inbounds X_pölli[i,j],Y_pölli[i,j] = func_cartesis(R_sektori[i,j],Θ_sektori_uus[i])
                end
            end
            #println("X_pölli[1,1] = ", X_pölli[1,1])
            #println("Y_pölli[1,1] = ", Y_pölli[1,1])
            # Lasketaan terän uusi asema
            @fastmath tera_asema_uus = tera_asema - (t0 / (2*pi) * ω*ii)
            # Muutetaan uus asema karteesiseen koordinaatistoon
            @simd for j in range(1,length(tera_asema))
                        @fastmath @inbounds X_terä[j,1],Y_terä[j,1] = func_cartesis(tera_asema_uus[j,1],tera_kulma)
                end
            # Poistetaan edellinen kuvaaja
            profile[:remove]()
            tera[:remove]()
            # Asetetaan datapisteille uudet asemat
            if k == 1
                profile = ax[:plot_wireframe](X_pölli, Y_pölli, Z_pölli)
            else
                profile = ax[:plot_surface](X_pölli, Y_pölli, Z_pölli)
            end
            tera = ax[:plot_wireframe](X_terä,Y_terä, Z_terä, color="r")

        else
            # Katotaan kuinka moni piste menee terälinjan ohi per askel
            @fastmath Θ_sektori_edellinen = Θ_sektori_uus.-ω
            #println("ii = ", ii)
            #println("Θ_sektori_uus[1,1] = ", Θ_sektori_uus[1,1])
            #println("Θ_sektori_edellinen[1,1] = ", Θ_sektori_edellinen[1,1])

            @simd for j in range(1,k)
                #println("j= ", j)
                piste = Int64[]
                Fc_sorvaus = Float64[]
                Fv_sorvaus = Float64[]

                @simd for i in range(1,n+1)
                    @fastmath @inbounds a1 = func_cartesis(R_sektori[i,j],Θ_sektori_edellinen[i]-tera_kulma)
                    @fastmath @inbounds a2 = func_cartesis(R_sektori[i,j],Θ_sektori_uus[i]-tera_kulma)
                    if sign(a1[1]) == sign(a2[1]) == 1 && sign(a1[2]) != sign(a2[2])
                        push!(piste,i)
                    end
                end
                #println("piste = ", piste)
                #println(length(piste))
                if length(piste) == 0 #Jos yksikään piste ei mene terälinjan ohi, skipataan voimien laskenta

                else
                    # Ehto joka tsekkaa onko terä pöllin sisäpuolella IF TRUE -> lasketaan leikkausvoima ja lasketaan pöllin pisteelle uusi koordinaatti
                    @simd for l in range(1,length(piste))
                        if tera_asema_uus[j] <= R_sektori[piste[l],j]
                            if j == k #Skipataan pöllin viimeisen pisteen kohdalla voiman laskenta, ettei lasketa ylimääräistä voimaa
                                # Lasketaan leikkauspaksuus
                                @fastmath @inbounds t0_temp = R_sektori[piste[l],j] - tera_asema_uus[j]
                                # Lasketaan leikatun pisteen uusi asema
                                @fastmath @inbounds R_sektori[piste[l],j] = R_sektori[piste[l],j] - t0_temp
                            else
                                # Lasketaan leikkauspaksuus
                                local t0_temp::Float64
                                @fastmath @inbounds t0_temp = R_sektori[piste[l],j] - tera_asema_uus[j]
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
                                local Fc_sorvaus_i::Float64 #Leikkauspinnan suuntainen voimakomponentti
                                local Fv_sorvaus_i::Float64 #Leikkauspintaan kohtisuoraan oleva voimakomponentti

                                @fastmath Fc_sorvaus_i = func_Fc(β_friction,φ_1,α_rake,τ,w_terä,γ,t0_temp,R)
                                push!(Fc_sorvaus,Fc_sorvaus_i)
                                @fastmath Fv_sorvaus_i = func_Fv(Fc_sorvaus_i,β_friction,α_rake)
                                push!(Fv_sorvaus,Fv_sorvaus_i)

                                # Lasketaan leikatun pisteen uusi asema
                                @fastmath @inbounds R_sektori[piste[l],j] = R_sektori[piste[l],j] - t0_temp
                            end
                        else
                            Fc_sorvaus_i = 0.0
                            push!(Fc_sorvaus,Fc_sorvaus_i)
                            Fv_sorvaus_i = 0.0
                            push!(Fv_sorvaus,Fv_sorvaus_i)
                        end
                    end
                    # Tallennetaan teräpisteen voimatiedot
                    if j == k
                        #Skipataan viimeinen siivu, koska ei lasketa sen voimia
                    else
                        #=if ii == 3
                            Fc_sorvaus_kootut[j] = copy(Fc_sorvaus)
                            Fv_sorvaus_kootut[j] = copy(Fv_sorvaus)
                        else=#
                            append!(Fc_sorvaus_kootut[j],Fc_sorvaus)
                            append!(Fv_sorvaus_kootut[j],Fv_sorvaus)
                        #end
                    end
                end
            end

            # Lasketaan pöllin uusi kulma-asema
            func_Θ_sektori_uus_serial!(Θ_sektori_uus,Θ_sektori,ω,ii)
            # Muutetaan uus asema karteesiseen koordinaatistoon
            @simd for j in range(1,k)
                @simd for i in 1:(n+1)
                    @fastmath @inbounds X_pölli[i,j],Y_pölli[i,j] = func_cartesis(R_sektori[i,j],Θ_sektori_uus[i])
                end
            end
            #println("X_pölli[1,1] = ", X_pölli[1,1])
            #println("Y_pölli[1,1] = ", Y_pölli[1,1])
            # Lasketaan terän uusi asema
            @fastmath tera_asema_uus = tera_asema - (t0 / (2*pi) * ω*ii)
            # Muutetaan uus asema karteesiseen koordinaatistoon
            @simd for j in range(1,length(tera_asema))
                @fastmath @inbounds X_terä[j,1],Y_terä[j,1] = func_cartesis(tera_asema_uus[j,1],tera_kulma)
            end
            ##################################
            # Pöllin dynamiikan laskenta
            ##################################
            # Laskee kaikkien profiilien pinta-alat
            A = func_area_elements(R_sektori, Θ_sektori, n, k)
            # Laskee kaikkien elementtien tilavuudet ja tallentaa ne matriisiin
            volume_elements = func_volume_elements(Z_pölli, A, n, k)
            # Laskee koko pöllin tilavuuden
            volume_total = func_volume_total(volume_elements)
            # laskee kaikkien elementtien tilavuuskeskiöt ja tallentaa ne Array of Arrays:iin
            CoV_elements = func_CoV_elements(X_pölli, Y_pölli, Z_pölli, A, n, k)
            # Laskee koko pöllin tilavuuskeskiön
            CoV_pölli = func_CoV_pölli(CoV_elements, volume_elements, volume_total, n, k)
            # Laskee kaikkien elementtien inertiat
            I_elements = func_inertias_elements(CoV_elements, X_pölli, Y_pölli, Z_pölli, ρ, n, k)
            # Laskee koko pöllin inertian
            I_total = func_inertias_pölli(CoV_elements, I_elements, volume_elements, ρ, CoV_pölli[1], CoV_pölli[2], CoV_pölli[3], n, k)
            ##################################
            # Kuvaajan piirto
            ##################################
            # Poistetaan edellinen kuvaaja
            profile[:remove]()
            tera[:remove]()
            # Asetetaan datapisteille uudet asemat
            if k == 1
                profile = ax[:plot_wireframe](X_pölli, Y_pölli, Z_pölli)
            else
                profile = ax[:plot_surface](X_pölli, Y_pölli, Z_pölli)
            end
            tera = ax[:plot_wireframe](X_terä,Y_terä, Z_terä, color="r")
        end
        ii = ii + 1
        fig[:canvas][:update]()
        fig[:canvas][:flush_events]()
        #sleep(0.001) #Julian oma sleep komento. Minimiaika on 1 ms
    end
    #return Fc_sorvaus_kootut, Fv_sorvaus_kootut
    return nothing
end #Function end

#Fc_data,Fv_data = sorvaus(R_alku,R_purilas,t0)
sorvaus(R_alku,R_purilas,t0,R_sektori)
#EOF
