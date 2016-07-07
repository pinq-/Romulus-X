# Funktiot joita käytetään leikkausvoimien laskemisessa 3D-tapauksessa
######################################################################
# Laskee pöllin profiilien z-koordinaatit. Serial == yhdellä prosessorilla.
function func_z_serial(k::Int64,w::Float64)
    wi::Float64 = w / (k-1)
    z = Array(Float64,k)
    @simd for i in range(1,k)
        if i == 1
            z[i] = 0.0
        else
            @fastmath z[i] = wi*(i-1)
        end
    end
    return z
end

#######################################
# Laskee pöllin R_sektorin & Θ_sektorin. Serial == yhdellä prosessorilla
function func_sektorit_serial(n::Int64,k::Int64,R_alku::Float64)
    R_sektori = Array(Float64, n, k)
    Θ_sektori = Array(Float64, n)
    @simd for j in range(1,k)
        # Yhden sektorin kulma
        #Θ_sektori = 2*pi / n
        # Jaetaan profiili n:ään pisteeseen
        @simd for i in range(1,n)
            @fastmath @inbounds R_sektori[i,j] = R_alku * (rand(89:105)/100)
            @fastmath @inbounds Θ_sektori[i] = (2*pi / n) * i
        end
    end
    # Lisätään profiilien datan alkupisteitä, jotta plottauksessa profiili ei ole aukinainen
    unshift!(Θ_sektori,0.0)
    #Θ_sektori = vcat(zeros(1), Θ_sektori) #Käytä vcatia jos on tarvetta määritellä profiilien pisteille eri kulmat. Nyt kaikki profiilit jaettu samoihin kulmiin
    R_sektori = vcat(R_sektori[n,:], R_sektori)

    return Θ_sektori, R_sektori
end

###################################
# Laskee pöllin pisteiden karteesiset koordinaatit. Serial = yhdellä prosessorilla
function func_XY_pölli_serial!(X,Y,R_sektori, Θ_sektori)
    @simd for j in range(1,k)
        @simd for i in 1:(n+1)
            @fastmath @inbounds X[i,j],Y[i,j] = func_cartesis(R_sektori[i,j],Θ_sektori[i])
        end
    end
    return X,Y
end

####################################
# Laskee pöllin uuden kulma-aseman
function func_Θ_sektori_uus_serial!(Θ_sektori_uus,Θ_sektori,ω,ii)
    @simd for j in range(1,size(Θ_sektori,2))
        @simd for i in range(1,size(Θ_sektori,1))
            @fastmath @inbounds Θ_sektori_uus[i,j] = Θ_sektori[i,j]-(ω*ii)
        end
    end
    return Θ_sektori_uus
end

####################################
# Laskee terän uuden aseman
function func_tera_asema_uus_serial!(tera_asema_uus::Array{Float64},tera_asema::Array{Float64},ω::Float64,ii::Int64)
    @simd for j in range(1,size(tera_asema,2))
        @simd for i in range(1,size(tera_asema,1))
            @fastmath @inbounds tera_asema_uus[i,j] = tera_asema[i,j]-(t0 / (2*pi) * ω*ii)
        end
    end
    return tera_asema_uus
end


#####################################################################
# PARALLEL PROCESSING FUNCTIONS
#####################################################################
# Laskee pöllin profiilien pisteet usealla prosessorilla
function func_R_sektori_parallel!(R_sektori)
    @sync begin
        for p in procs(R_sektori)
            @async remotecall_wait(R_sektori_shared_chunk!, p, R_sektori)
        end
    end
    #R_sektori
    return nothing
end

# Here's a convenience wrapper for a SharedArray implementation
@everywhere R_sektori_shared_chunk!(q) = R_sektori_chunk!(q, myrange(q)...)

@everywhere function R_sektori_chunk!(q, irange, jrange)
    #@show (irange, jrange)  # display so we can see what's happening
    @simd for j in jrange
        @simd for i in irange
                @fastmath @inbounds q[i,j] = R_alku * (rand(89:105)/100)
        end
    end
    #q
    return nothing
end

####################################
# Laskee Θ_sektorin usealla prosessorilla
function func_Θ_sektori_parallel!(Θ_sektori)
    @sync begin
        for p in procs(Θ_sektori)
            @async remotecall_wait(Θ_sektori_shared_chunk!, p, Θ_sektori)
        end
    end
    #Θ_sektori
    return nothing
end
# Here's a convenience wrapper for a SharedArray implementation
@everywhere Θ_sektori_shared_chunk!(q) = Θ_sektori_chunk!(q, myrangerows(q)...)

@everywhere function Θ_sektori_chunk!(q, irange, jrange)
    #@show (irange, jrange)  # display so we can see what's happening
    @simd for j in jrange
        @simd for i in irange
            if i == 1
                q[i,j] = 0.0
            else
                @fastmath @inbounds q[i,j] = (2*pi / n) * (i-1)
            end
        end
    end
    #q
    return nothing
end

######################################
# Laskee profiilien pisteiden karteesiset koordinaatit usealla prosessorilla
@everywhere function func_XY_pölli_parallel!(X,Y, R_sektori, Θ_sektori)
    @sync begin
        for p in procs(X)
            @async remotecall_wait(XY_shared_chunk!, p, X, Y, R_sektori, Θ_sektori)
        end
    end
    #X,Y
    return nothing
end

# Here's a convenience wrapper for a SharedArray implementation
@everywhere XY_shared_chunk!(q,u, R_sektori, Θ_sektori) = XY_chunk!(q, u, R_sektori, Θ_sektori, myrange(q)...)

@everywhere function XY_chunk!(q, u, R_sektori, Θ_sektori, irange, jrange)
    #@show (irange, jrange)  # display so we can see what's happening
    @simd for j in jrange
        @simd for i in irange
            @fastmath @inbounds q[i,j],u[i,j] = func_cartesis(R_sektori[i,j],Θ_sektori[i])
        end
    end
    #q, u
    return nothing
end

#########################################
# Laskee profiilien uusien pisteiden karteesiset koordinaatit usealla prosessorilla
# Käytetään while loopin sisällä
@everywhere function func_XY_pölli_uudet_parallel!(X::SharedArray{Float64},Y::SharedArray{Float64}, R_sektori::SharedArray{Float64}, Θ_sektori_uus::SharedArray{Float64})
    @sync begin
        for p in procs(X)
            @async remotecall_wait(XY_uudet_shared_chunk!, p, X, Y, R_sektori, Θ_sektori_uus)
        end
    end
    #X,Y
    return nothing
end

# Here's a convenience wrapper for a SharedArray implementation
@everywhere XY_uudet_shared_chunk!(q::SharedArray{Float64},u::SharedArray{Float64}, R_sektori::SharedArray{Float64}, Θ_sektori_uus::SharedArray{Float64}) = XY_uudet_chunk!(q, u, R_sektori, Θ_sektori_uus, myrange(q)...)

@everywhere function XY_uudet_chunk!(q::SharedArray{Float64}, u::SharedArray{Float64}, R_sektori::SharedArray{Float64}, Θ_sektori_uus::SharedArray{Float64}, irange::UnitRange{Int64}, jrange::UnitRange{Int64})
    #@show (irange, jrange)  # display so we can see what's happening
    @simd for j in jrange
        @simd for i in irange
            @fastmath @inbounds q[i,j],u[i,j] = func_cartesis(R_sektori[i,j],Θ_sektori_uus[i])
        end
    end
    #q, u
    return nothing
end


@everywhere function func_XY_pölli_uudet_parallel2!(X_pölli::SharedArray{Float64},Y_pölli::SharedArray{Float64}, R_sektori::SharedArray{Float64}, Θ_sektori_uus::SharedArray{Float64})
    @sync @parallel for j in range(1,k)
        @simd for i in range(1,n+1)
            @fastmath @inbounds X_pölli[i,j],Y_pölli[i,j] = func_cartesis(R_sektori[i,j],Θ_sektori_uus[i])
        end
    end
    return nothing
end
################################
# Terän pituuden laskenta - parallel
@everywhere function func_z_terä_parallel!(z_terä, tera_asema, w_terä, tera_apukulma)
    @sync begin
        for p in procs(z_terä)
            @async remotecall_wait(z_terä_shared_chunk!, p, z_terä, tera_asema, w_terä, tera_apukulma)
        end
    end
    z_terä
end

# Here's a convenience wrapper for a SharedArray implementation
@everywhere z_terä_shared_chunk!(z_terä, tera_asema, w_terä, tera_apukulma) = z_terä_chunk!(z_terä, tera_asema, w_terä, tera_apukulma, myrangerows(z_terä)...)

@everywhere function z_terä_chunk!(z_terä, tera_asema, w_terä, tera_apukulma, irange::UnitRange{Int64}, jrange::UnitRange{Int64})
    #@show (irange, jrange)  # display so we can see what's happening
    @simd for j in jrange
        @simd for i in irange
            if i == 1
                z_terä[i,j] = 0.0
            elseif i < length(tera_asema)
                @fastmath @inbounds z_terä[i,j] = w_terä*(i-1)
                @fastmath @inbounds tera_asema[i,j] = tera_asema[1,j] + (z_terä[i,j] * tan(tera_apukulma))
            else
                @fastmath @inbounds z_terä[i,j] = w_terä*(i-1)
            end
        end
    end
    z_terä
end

####################################
# Laskee Θ_sektori_uus usealla prosessorilla
function func_Θ_sektori_uus_parallel!(Θ_sektori_uus, Θ_sektori, ω, ii)
    @sync begin
        for p in procs(Θ_sektori_uus)
            @async remotecall_wait(Θ_sektori_uus_shared_chunk!, p, Θ_sektori_uus, Θ_sektori, ω, ii)
        end
    end
    #Θ_sektori_uus
    return nothing
end
# Here's a convenience wrapper for a SharedArray implementation
@everywhere Θ_sektori_uus_shared_chunk!(q, Θ_sektori, ω, ii) = Θ_sektori_uus_chunk!(q, Θ_sektori, ω, ii, myrangerows(q)...)

@everywhere function Θ_sektori_uus_chunk!(q, Θ_sektori, ω, ii, irange, jrange)
    #@show (irange, jrange)  # display so we can see what's happening
    @simd for j in jrange
        @simd for i in irange
            @fastmath @inbounds q[i,j] = Θ_sektori[i,j]-(ω*ii)
        end
    end
    #q
    return nothing
end
