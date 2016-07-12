# Lisätään Julialle koneen ytimet käyttöön
if nprocs() == 1
    addprocs(CPU_CORES - 1)
elseif nprocs() == CPU_CORES

else
    rmprocs(workers())
    addprocs(CPU_CORES - 1)
end


@everywhere include("leikkausvoimien_laskenta.jl")
include("3D_laskentafunktiot.jl")
@everywhere include("lahto_arvot.jl")
@everywhere include("apufunktiot.jl")
include("parallel_tutorial.jl")

using PyCall
using PyPlot
@pyimport matplotlib.cm as cm

# Jaetaan pölli k:aan profiiliin
@everywhere const k = 486
# Pisteiden lukumäärä
@everywhere const n = 359
# Kulma, jonka verran pöllin profiili pyörii per askel
@everywhere const ω = func_radians(60.0) #[rad]


#Animaatiofunktio
function sorvaus(R_alku::Float64, R_purilas::Float64, t0::Float64)

    R_sektori::SharedArray{Float64} = SharedArray(Float64,n+1,k, init = R_sektori -> R_sektori[Base.localindexes(R_sektori)] = myid())
    Θ_sektori::SharedArray{Float64} = SharedArray(Float64, n+1, init = Θ_sektori -> Θ_sektori[Base.localindexes(Θ_sektori)] = myid())
    X_pölli::SharedArray{Float64} = SharedArray(Float64, n+1, k, init = X_pölli -> X_pölli[Base.localindexes(X_pölli)] = myid())
    Y_pölli::SharedArray{Float64} = SharedArray(Float64, n+1, k, init = Y_pölli -> Y_pölli[Base.localindexes(Y_pölli)] = myid())
    Θ_sektori_uus::SharedArray{Float64} = SharedArray(Float64, n+1, init = Θ_sektori_uus -> Θ_sektori_uus[Base.localindexes(Θ_sektori_uus)] = myid())
    Θ_sektori_edellinen::SharedArray{Float64} = SharedArray(Float64, n+1, init = Θ_sektori_edellinen -> Θ_sektori_edellinen[Base.localindexes(Θ_sektori_edellinen)] = myid())
    Fc_sorvaus_kootut::SharedArray{Float64} = SharedArray(Float64,1,k-1)

    # Lasketaan profiilien z-koordinaatit
    z = func_z_serial(k,w)

    # Kuvataan profiili polaarikoordinaateissa. Jokainen piste on muotoa (Θ,r)
    # Jaetaan profiili n+1:een pisteeseen
    func_R_sektori_parallel!(R_sektori)
    func_Θ_sektori_parallel!(Θ_sektori)
    # Lisätään profiilien dataan alkupisteitä, jotta plottauksessa profiili ei ole aukinainen
    for j in range(1,k)
        R_sektori[1,j] = R_sektori[n+1,j]
    end


    # Muutetaan profiilien pisteet karteesiseen koordinaatistoon
    func_XY_pölli_parallel!(X_pölli,Y_pölli, R_sektori, Θ_sektori)

#=    # Määrittää kuvaajan
    fig = plt[:figure]()
    ax = fig[:add_subplot](111, projection="3d")
    # Määritellään koordinaatiston rajat
    ax[:set_xlim3d](-R_alku-0.53,R_alku+0.53)
    ax[:set_zlim3d](0,w+0.1)
    ax[:set_ylim3d](-R_alku-0.53,R_alku+0.53)
=#

    # Terän paikka koordinaatistossa
    tera_asema = Array(Float64, k)
    tera_asema_uus::SharedArray{Float64} = SharedArray(Float64, size(tera_asema,1))
    z_terä = Array(Float64, length(tera_asema))
    #z_terä = SharedArray(Float64, length(tera_asema),init = X_pölli -> X_pölli[Base.localindexes(X_pölli)] = myid())
    w_terä = w / (length(tera_asema)-1)
    #Terän päätepisteet. Jos eri niin looppi laskee terän pisteille etäisyyden origosta
    tera_asema[1] = 0.13
    tera_asema[k] = 0.13
    tera_apukulma::Float64 = atan( (tera_asema[k] - tera_asema[1])/w )
    # Lasketaan terän pituus
    #func_z_terä_parallel!(z_terä, tera_asema, w_terä, tera_apukulma)
    @simd for i in range(1,length(tera_asema))
        if i == 1
            z_terä[i] = 0.0
        elseif i < length(tera_asema)
            @fastmath @inbounds z_terä[i] = w_terä*(i-1)
            @fastmath @inbounds tera_asema[i] = tera_asema[1] + (z_terä[i] * tan(tera_apukulma))
        else
            @fastmath @inbounds z_terä[i] = w_terä*(i-1)
        end
    end
    #tera_kulma::Float64 = func_radians(0) #[rad]
    @everywhere tera_kulma = func_radians(0.0) #[rad]

    # Muutetaan terän pisteet karteesiseen koordinaatistoon
    X_terä = Array(Float64, length(tera_asema))
    Y_terä = Array(Float64, length(tera_asema))
    @simd for j in range(1,length(tera_asema))
            @fastmath @inbounds X_terä[j,1],Y_terä[j,1] = func_cartesis(tera_asema[j,1],tera_kulma)
    end



    # Plotataan profiili, lasketaan joka aika-askeleella pisteille uusi kulma-asema ja päivitetään plottaus
    ii::Int64 = 1
    while R_sektori[1,1] >= R_purilas
        if ii == 1
#=            ax[:view_init](elev=90, azim=-90)
            fig[:canvas][:draw]()
            # Piirretään pöllin alkuasema
            if k == 1
                profile = ax[:plot_wireframe](X_pölli, Y_pölli, z)
            else
                profile = ax[:plot_surface](X_pölli, Y_pölli, z)
            end
            # Terän alkuasema
            tera = ax[:plot_wireframe](X_terä,Y_terä, z_terä, color="r")
=#
        elseif ii == 2
            # Lasketaan pöllin uusi kulma-asema
            #func_Θ_sektori_uus_serial!(Θ_sektori_uus,Θ_sektori,ω,ii)
            func_Θ_sektori_uus_parallel!(Θ_sektori_uus, Θ_sektori, ω, ii)

            # Lasketaan terän uusi asema
            @simd for j in range(1,size(tera_asema,2))
                @simd for i in range(1,size(tera_asema,1))
                    @fastmath @inbounds tera_asema_uus[i,j] = tera_asema[i,j]-(t0 / (2*pi) * ω*ii)
                end
            end


            # Muutetaan pöllin asema karteesiseen koordinaatistoon
            func_XY_pölli_uudet_parallel!(X_pölli,Y_pölli, R_sektori, Θ_sektori_uus)
            #func_XY_pölli_serial!(X_pölli,Y_pölli,R_sektori, Θ_sektori_uus)

            # Muutetaan terän uus asema karteesiseen koordinaatistoon
            @simd for j in range(1,length(tera_asema))
                        @fastmath @inbounds X_terä[j,1],Y_terä[j,1] = func_cartesis(tera_asema_uus[j,1],tera_kulma)
                end
#=           # Poistetaan edellinen kuvaaja
            profile[:remove]()
            tera[:remove]()
            # Asetetaan datapisteille uudet asemat
            if k == 1
                profile = ax[:plot_wireframe](X_pölli, Y_pölli, z)
            else
                profile = ax[:plot_surface](X_pölli, Y_pölli, z)
            end
            tera = ax[:plot_wireframe](X_terä,Y_terä, z_terä, color="r")
=#
        else
            # Katotaan kuinka moni piste menee terälinjan ohi per askel
            for i in range(1,n+1)
                @fastmath @inbounds Θ_sektori_edellinen[i,1] = Θ_sektori_uus[i,1]-ω
            end

            piste = Int64[]
            @simd for i in range(1,n+1)
                @fastmath @inbounds a1 = func_cartesis(R_sektori[i,1],Θ_sektori_edellinen[i]-tera_kulma)
                @fastmath @inbounds a2 = func_cartesis(R_sektori[i,1],Θ_sektori_uus[i]-tera_kulma)
                if sign(a1[1]) == sign(a2[1]) == 1 && sign(a1[2]) != sign(a2[2])
                    push!(piste,i)
                end
            end
            piste_temp::SharedArray{Int64}
            piste_temp = SharedArray(Int64, length(piste),k)
            Fc_temp::SharedArray{Float64}
            Fc_temp = SharedArray(Float64, length(piste),k-1)

            func_piste_parallel!(piste_temp, R_sektori, Θ_sektori_edellinen, Θ_sektori_uus, tera_kulma)
            func_Fc_parallel!(piste_temp, Fc_temp, Fc_sorvaus_kootut, tera_asema_uus, R_sektori, w_terä, ii)
            Fc_sorvaus_kootut = vcat(Fc_sorvaus_kootut,Fc_temp)

#=            @simd for j in range(1,k)
                piste = Int64[]
                Fc_sorvaus = Float64[]

                @simd for i in range(1,n+1)
                    @fastmath @inbounds a1 = func_cartesis(R_sektori[i,j],Θ_sektori_edellinen[i]-tera_kulma)
                    @fastmath @inbounds a2 = func_cartesis(R_sektori[i,j],Θ_sektori_uus[i]-tera_kulma)
                    if sign(a1[1]) == sign(a2[1]) == 1 && sign(a1[2]) != sign(a2[2])
                        push!(piste,i)
                    end
                end
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
                                local Fc_sorvaus_i::Float64
                                @fastmath Fc_sorvaus_i = func_Fc(β_friction,φ_1,α_rake,τ,w_terä,γ,t0_temp,R)
                                push!(Fc_sorvaus,Fc_sorvaus_i)

                                # Lasketaan leikatun pisteen uusi asema
                                @fastmath @inbounds R_sektori[piste[l],j] = R_sektori[piste[l],j] - t0_temp
                            end
                        else
                            Fc_sorvaus_i = 0.0
                            push!(Fc_sorvaus,Fc_sorvaus_i)
                        end
                    end
                    # Tallennetaan teräpisteen voimatiedot
                    if j == k
                        #Skipataan viimeinen siivu, koska ei lasketa sen voimia
                    else
                        if ii == 3
                            Fc_sorvaus_kootut[j] = copy(Fc_sorvaus)
                        else
                            append!(Fc_sorvaus_kootut[j],Fc_sorvaus)
                        end
                    end
                end
            end=#

            # Lasketaan pöllin uusi kulma-asema
            #func_Θ_sektori_uus_serial!(Θ_sektori_uus,Θ_sektori,ω,ii)
            func_Θ_sektori_uus_parallel!(Θ_sektori_uus, Θ_sektori, ω, ii)

            # Lasketaan terän uusi asema
            #@fastmath tera_asema_uus = tera_asema - (t0 / (2*pi) * ω*ii)
            @simd for j in range(1,size(tera_asema,2))
                @simd for i in range(1,size(tera_asema,1))
                    @fastmath @inbounds tera_asema_uus[i,j] = tera_asema[i,j]-(t0 / (2*pi) * ω*ii)
                end
            end

            # Muutetaan pöllin uus asema karteesiseen koordinaatistoon
            func_XY_pölli_uudet_parallel!(X_pölli,Y_pölli, R_sektori, Θ_sektori_uus)
            #func_XY_pölli_serial!(X_pölli,Y_pölli,R_sektori, Θ_sektori_uus)

            # Muutetaan terän uus asema karteesiseen koordinaatistoon
            @simd for j in range(1,length(tera_asema))
                        @fastmath @inbounds X_terä[j,1],Y_terä[j,1] = func_cartesis(tera_asema_uus[j,1],tera_kulma)
                end
#=            # Poistetaan edellinen kuvaaja
            profile[:remove]()
            tera[:remove]()
            # Asetetaan datapisteille uudet asemat
            if k == 1
                profile = ax[:plot_wireframe](X_pölli, Y_pölli, z)
            else
                profile = ax[:plot_surface](X_pölli, Y_pölli, z)
            end
            tera = ax[:plot_wireframe](X_terä,Y_terä, z_terä, color="r")
=#
        end
        ii = ii + 1
        #fig[:canvas][:update]()
        #fig[:canvas][:flush_events]()
        #sleep(0.001) #Julian oma sleep komento. Minimiaika on 1 ms

    end
    return Fc_sorvaus_kootut
end #Function end

#data = sorvaus(R_alku, R_purilas, t0, R_sektori, theta_sektori, theta_sektori_uus, X_pölli, Y_pölli, Θ_sektori_edellinen)
#EOF
