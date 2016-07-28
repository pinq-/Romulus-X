# Pöllin dynamiikan laskentafunktiot #
######################################

# Laskee R_sektori[i,j] ja R_sektori[i+1,j] välisten pisteiden etäisyyden
function func_d(R1::Float64, R2::Float64, Θ1::Float64, Θ2::Float64)
    d = sqrt( R1^2 + R2^2 - 2*R1*R2*cos(Θ2 - Θ1) ) #[m]
end

# Heron's yhtälö kolmion pinta-alalle
function func_triangle_area(a::Float64, b::Float64, c::Float64)
    s = (a+b+c) / 2
    area = sqrt( s*(s-a)*(s-b)*(s-c) ) #[m^2]
end

# Laskee yhden profiilin kokonaispinta-alan
function func_profile_area(R_sektori::Array{Float64},Θ_sektori::Array{Float64})
        area_i = Float64[]
        for i in range(1,size(R_sektori,1)-1)
                d = func_d(R_sektori[i], R_sektori[i+1], Θ_sektori[i], Θ_sektori[i+1])
                push!(area_i,func_triangle_area(R_sektori[i], R_sektori[i+1], d))
        end
        area_profile = sum(area_i,1) #[m^2]
    return area_profile
end

# Laskee yhden elementin tilavuuden
# Truncated triangle pyramid volume equation
function func_volume_i(h::Float64,A1::Float64,A2::Float64)
    volume_i = h/3 * (A1 + sqrt( A1 * A2 ) + A2) #[m^3]
end

# Laskee kaikkien elementtien tilavuudet ja tallentaa ne matriisiin
function func_volume_elements(Z_pölli::Array{Float64}, A::Array{Float64}, n::Int64, k::Int64)
    volume_elements = Array(Float64, n, k-1)

    @simd for j in range(1,k-1)
        h::Float64 = Z_pölli[j+1] - Z_pölli[j]
        @simd for i in range(1,n)
            @inbounds @fastmath volume_elements[i,j] = func_volume_i(h,A[i,j],A[i,j+1])
        end
    end
    return volume_elements
end

# Laskee kaikkien profiilien osien pinta-alat
function func_area_elements(R_sektori::Array{Float64}, Θ_sektori::Array{Float64}, n::Int64, k::Int64)
    c = Array(Float64, n, k)
    A = Array(Float64, n, k)
    @simd for j in range(1,k)
        @simd for i in range(1,n)
            @inbounds @fastmath c[i,j] = func_d(R_sektori[i,j], R_sektori[i+1,j], Θ_sektori[i], Θ_sektori[i+1])
            @inbounds @fastmath A[i,j] = func_triangle_area(R_sektori[i,j], R_sektori[i+1,j], c[i,j])
        end
    end
    return A
end

# Laskee koko pöllin tilavuuden
function func_volume_total(volume_elements::Array{Float64})
    # Lasketaan koko pöllin tilavuus
    volume_total::Float64 = sum(volume_elements) #[m^3]
end

# Laskee pöllin kokonaismassan
function func_pölli_massa(ρ::Float64, volume_total::Float64)
    massa = ρ * volume_total #[kg]
end

#########################################
# Elementin massakeskipisteen laskentaa #
#########################################
# Laskee yhden elementin painopisteen Z-koordinaatin
function func_CoV_Z_i(h::Float64, A1::Float64, A2::Float64)
    Z_i = (h/4) * (A1 + 2*sqrt(A1*A2)+3*A2) / (A1+sqrt(A1*A2)+A2)
end

# Laskee kolmion pintakeskiön X- ja Y-koordinaatit
function func_triangle_centroid(X_pölli::Array{Float64}, Y_pölli::Array{Float64}, i::Int64, j::Int64)
    # Kolmannen pisteen koordinatit
    Ox::Float64 = 0.0
    Oy::Float64 = 0.0
    # pintakeskiön laskenta
    centroid_x = (Ox + X_pölli[i,j] + X_pölli[i+1,j]) / 3.0
    centroid_y = (Oy + Y_pölli[i,j] + Y_pölli[i+1,j]) / 3.0

    return centroid_x,centroid_y
end

function func_triangle_centroid2(X_pölli::Float64, Y_pölli::Float64, X_pölli2::Float64, Y_pölli2::Float64)
    # Kolmannen pisteen koordinatit
    Ox::Float64 = 0.0
    Oy::Float64 = 0.0
    # pintakeskiön laskenta
    centroid_x = (Ox + X_pölli + X_pölli2) / 3.0
    centroid_y = (Oy + Y_pölli + Y_pölli2) / 3.0

    return centroid_x,centroid_y
end

# Laskee yhden elementin painopisteen koordinaatit
function func_CoV_element(X_pölli::Array{Float64}, Y_pölli::Array{Float64}, Z_pölli::Array{Float64}, i::Int64, j::Int64, A::Array{Float64})
    p1 = Array(Float64,3,1)
    p2 = Array(Float64,3,1)
    vektori = Array(Float64,3,1)
    cvektori = Array(Float64,3,1)
    CoV_element = Array(Float64,3,1)
    # Ratkaistaan kolmioiden pintakeskiöt
    p1[1,1],p1[2,1] = func_triangle_centroid2(X_pölli[i,j], Y_pölli[i,j], X_pölli[i+1,j], Y_pölli[i+1,j])
#    p1[1,1],p1[2,1] = func_triangle_centroid(X_pölli, Y_pölli, i, j)
    p1[3,1] = Z_pölli[j]
    p2[1,1],p2[2,1] = func_triangle_centroid2(X_pölli[i,j+1], Y_pölli[i,j+1], X_pölli[i+1,j+1], Y_pölli[i+1,j+1])
#    p2[1,1],p2[2,1] = func_triangle_centroid(X_pölli, Y_pölli, i, j+1)
    p2[3,1] = Z_pölli[j+1]

    h::Float64 = p2[3,1] - p1[3,1]
    # Laskee elementin painopisteen Z-koordinaatin
    Z_i = func_CoV_Z_i(h, A[i,j], A[i,j+1])

    # Vektori, joka kulkee kolmioiden pintakeskiöiden kautta.
    #vektori = p2 - p1
    for l in range(1,3)
        vektori[l] = p2[l] - p1[l]
    end
    # Vektori, joka on samansuuntainen edellisen vektorin kanssa. Kulkee painopisteeseen.
    #cvektori = vektori * (Z_i/h)
    for l in range(1,3)
        cvektori[l] = vektori[l] * (Z_i/h)
    end
    # Ratkaistaan ensimmäisen kolmion pintakeskiön ja cvektorin avulla painopisteen koordinaatit
    #CoV_element = cvektori + p1
    for l in range(1,3)
        CoV_element[l] = cvektori[l] + p1[l]
    end
return CoV_element
end

function func_CoV_element2(X_pölli11, Y_pölli11, X_pölli21, Y_pölli21, X_pölli12, Y_pölli12, X_pölli22, Y_pölli22, Z_pölli1, Z_pölli2, A11, A12)
    p1 = Array(Float64,3,1)
    p2 = Array(Float64,3,1)
    vektori = Array(Float64,3,1)
    cvektori = Array(Float64,3,1)
    CoV_element = Array(Float64,3,1)
    # Ratkaistaan kolmioiden pintakeskiöt
    p1[1,1],p1[2,1] = func_triangle_centroid2(X_pölli11, Y_pölli11, X_pölli21, Y_pölli21)
#    p1[1,1],p1[2,1] = func_triangle_centroid(X_pölli, Y_pölli, i, j)
    p1[3,1] = Z_pölli1
    p2[1,1],p2[2,1] = func_triangle_centroid2(X_pölli12, Y_pölli12, X_pölli22, Y_pölli22)
#    p2[1,1],p2[2,1] = func_triangle_centroid(X_pölli, Y_pölli, i, j+1)
    p2[3,1] = Z_pölli2

    h::Float64 = p2[3,1] - p1[3,1]
    # Laskee elementin painopisteen Z-koordinaatin
    Z_i = func_CoV_Z_i(h, A11, A12)

    # Vektori, joka kulkee kolmioiden pintakeskiöiden kautta.
    #vektori = p2 - p1
    for l in range(1,3)
        vektori[l] = p2[l] - p1[l]
    end
    # Vektori, joka on samansuuntainen edellisen vektorin kanssa. Kulkee painopisteeseen.
    #cvektori = vektori * (Z_i/h)
    for l in range(1,3)
        cvektori[l] = vektori[l] * (Z_i/h)
    end
    # Ratkaistaan ensimmäisen kolmion pintakeskiön ja cvektorin avulla painopisteen koordinaatit
    #CoV_element = cvektori + p1
    for l in range(1,3)
        CoV_element[l] = cvektori[l] + p1[l]
    end
return CoV_element
end


# laskee kaikkien elementtien tilavuuskeskiöt ja tallentaa ne Array of Arrays:iin
function func_CoV_elements(X_pölli::Array{Float64}, Y_pölli::Array{Float64}, Z_pölli::Array{Float64}, A::Array{Float64}, n::Int64, k::Int64)
    CoV_elements = Array(Array{Float64},n,k-1)

    @simd for j in range(1,k-1)
      @simd for i in range(1,n)
#        @inbounds @fastmath CoV_elements[i,j] = func_CoV_element(X_pölli, Y_pölli, Z_pölli, i, j, A)
        @inbounds @fastmath CoV_elements[i,j] = func_CoV_element2(X_pölli[i,j], Y_pölli[i,j], X_pölli[i+1,j], Y_pölli[i+1,j], X_pölli[i,j+1], Y_pölli[i,j+1], X_pölli[i+1,j+1], Y_pölli[i+1,j+1], Z_pölli[j], Z_pölli[j+1], A[i,j], A[i,j+1])
      end
    end

    return CoV_elements
end

# Laskee koko pöllin tilavuuskeskiön
function func_CoV_pölli(CoV::Array{Array{Float64}}, volume_elements::Array{Float64}, volume_total::Float64, n::Int64, k::Int64)
    xi_sum = 0.0
    yi_sum = 0.0
    zi_sum = 0.0
    @simd for j in range(1,k-1)
        @simd for i in range(1,n)
            @inbounds @fastmath xi_sum = CoV[i,j][1]*volume_elements[i,j] + xi_sum
            @inbounds @fastmath yi_sum = CoV[i,j][2]*volume_elements[i,j] + yi_sum
            @inbounds @fastmath zi_sum = CoV[i,j][3]*volume_elements[i,j] + zi_sum
        end
    end

    x_CoV = xi_sum / volume_total
    y_CoV = yi_sum / volume_total
    z_CoV = zi_sum / volume_total
    return x_CoV,y_CoV,z_CoV
end

########################
# Inertioiden laskenta #
########################

# Laskee yhden elementin inertiat elementin keskipisteen suhteen
function func_inertias_element(CoV::Array{Array{Float64}}, X_pölli::Array{Float64}, Y_pölli::Array{Float64}, Z_pölli::Array{Float64}, ρ::Float64, i::Int64, j::Int64)
    # Jaetaan elementti 8:aan tetrahedraan, lasketaan tetrahedrojen inertiat elementin keskipisteen suhteen ja summataan yhteen
    Dvec = Array(Array{Float64},8)
    Evec = Array(Array{Float64},8)
    Fvec = Array(Array{Float64},8)
    G = Array(Array{Float64},8)
    H = Array(Array{Float64},8)
    N =  Array(Array{Float64},8)
    volume_tetra = Array(Float64,8)
    Pxx =  Array(Float64,8)
    Pyy =  Array(Float64,8)
    Pzz =  Array(Float64,8)
    Pxy =  Array(Float64,8)
    Pxz =  Array(Float64,8)
    Pyz =  Array(Float64,8)

    origo_temp = CoV[i,j]
    # Jaetaan tetrahedroihin
    Dvec[1] = [0.0-origo_temp[1]; 0.0-origo_temp[2]; Z_pölli[j]-origo_temp[3]]
    Evec[1] = [X_pölli[i,j]-origo_temp[1]; Y_pölli[i,j]-origo_temp[2]; Z_pölli[j]-origo_temp[3]]
    Fvec[1] = [X_pölli[i+1,j]-origo_temp[1]; Y_pölli[i+1,j]-origo_temp[2]; Z_pölli[j]-origo_temp[3]]

    Dvec[2] = [0.0-origo_temp[1]; 0.0-origo_temp[2]; Z_pölli[j+1]-origo_temp[3]]
    Evec[2] = [X_pölli[i+1,j+1]-origo_temp[1]; Y_pölli[i+1,j+1]-origo_temp[2]; Z_pölli[j+1]-origo_temp[3]]
    Fvec[2] = [X_pölli[i,j+1]-origo_temp[1]; Y_pölli[i,j+1]-origo_temp[2]; Z_pölli[j+1]-origo_temp[3]]

    Dvec[3] = [0.0-origo_temp[1]; 0.0-origo_temp[2]; Z_pölli[j]-origo_temp[3]]
    Evec[3] = [X_pölli[i+1,j]-origo_temp[1]; Y_pölli[i+1,j]-origo_temp[2]; Z_pölli[j]-origo_temp[3]]
    Fvec[3] = [0.0-origo_temp[1]; 0.0-origo_temp[2]; Z_pölli[j+1]-origo_temp[3]]

    Dvec[4] = [X_pölli[i+1,j]-origo_temp[1]; Y_pölli[i+1,j]-origo_temp[2]; Z_pölli[j]-origo_temp[3]]
    Evec[4] = [X_pölli[i+1,j+1]-origo_temp[1]; Y_pölli[i+1,j+1]-origo_temp[2]; Z_pölli[j+1]-origo_temp[3]]
    Fvec[4] = [0.0-origo_temp[1]; 0.0-origo_temp[2]; Z_pölli[j+1]-origo_temp[3]]

    Dvec[5] = [X_pölli[i+1,j]-origo_temp[1]; Y_pölli[i+1,j]-origo_temp[2]; Z_pölli[j]-origo_temp[3]]
    Evec[5] = [X_pölli[i,j+1]-origo_temp[1]; Y_pölli[i,j+1]-origo_temp[2]; Z_pölli[j+1]-origo_temp[3]]
    Fvec[5] = [X_pölli[i+1,j+1]-origo_temp[1]; Y_pölli[i+1,j+1]-origo_temp[2]; Z_pölli[j+1]-origo_temp[3]]

    Dvec[6] = [X_pölli[i+1,j]-origo_temp[1]; Y_pölli[i+1,j]-origo_temp[2]; Z_pölli[j]-origo_temp[3]]
    Evec[6] = [X_pölli[i,j]-origo_temp[1]; Y_pölli[i,j]-origo_temp[2]; Z_pölli[j]-origo_temp[3]]
    Fvec[6] = [X_pölli[i,j+1]-origo_temp[1]; Y_pölli[i,j+1]-origo_temp[2]; Z_pölli[j+1]-origo_temp[3]]

    Dvec[7] = [0.0-origo_temp[1]; 0.0-origo_temp[2]; Z_pölli[j]-origo_temp[3]]
    Evec[7] = [X_pölli[i,j+1]-origo_temp[1]; Y_pölli[i,j+1]-origo_temp[2]; Z_pölli[j+1]-origo_temp[3]]
    Fvec[7] = [X_pölli[i,j]-origo_temp[1]; Y_pölli[i,j]-origo_temp[2]; Z_pölli[j]-origo_temp[3]]

    Dvec[8] = [0.0-origo_temp[1]; 0.0-origo_temp[2]; Z_pölli[j]-origo_temp[3]]
    Evec[8] = [0.0-origo_temp[1]; 0.0-origo_temp[2]; Z_pölli[j+1]-origo_temp[3]]
    Fvec[8] = [X_pölli[i,j+1]-origo_temp[1]; Y_pölli[i,j+1]-origo_temp[2]; Z_pölli[j+1]-origo_temp[3]]

    for l in range(1,8)
        # Lasketaan tetrahedran pohjan reunat
        G[l] = Evec[l] - Dvec[l]
        H[l] = Fvec[l] - Dvec[l]
        # Lasketaan tetrahedran pohjan normaali
        #N[l] = cross(G[l],H[l])
        N[l] = cross(H[l],G[l])
        # Lasketaan tetrahedran tilavuus
        volume_tetra[l] = dot(Dvec[l]/3,N[l]/2)
    end
    # Lasketaan elementin tilavuus (Pelkästään tsekkausta varten)
    #volume = sum(volume_tetra)
    #println(volume)

    # Lasketaan tetrahedrojen products of inertiat Pjk, missä j,k = x, y, z = 1, 2, 3
    for i in range(1,8)
        Pxx[i] = func_Pjk_tetra(ρ, volume_tetra, Dvec, Evec, Fvec, i, 1, 1)
        Pyy[i] = func_Pjk_tetra(ρ, volume_tetra, Dvec, Evec, Fvec, i, 2, 2)
        Pzz[i] = func_Pjk_tetra(ρ, volume_tetra, Dvec, Evec, Fvec, i, 3, 3)

        Pxy[i] = func_Pjk_tetra(ρ, volume_tetra, Dvec, Evec, Fvec, i, 1, 2)
        Pxz[i] = func_Pjk_tetra(ρ, volume_tetra, Dvec, Evec, Fvec, i, 1, 3)
        Pyz[i] = func_Pjk_tetra(ρ, volume_tetra, Dvec, Evec, Fvec, i, 2, 3)
    end

    # Lasketaan elementin inertiat
    Ixx = sum(Pyy) + sum(Pzz)
    Iyy = sum(Pxx) + sum(Pzz)
    Izz = sum(Pxx) + sum(Pyy)
    Iyz = -sum(Pyz)
    Ixz = -sum(Pxz)
    Ixy = -sum(Pxy)

    I_element = [Ixx Ixy Ixz; Ixy Iyy Iyz; Ixz Iyz Izz]

    return I_element
end

function func_inertias_element2(CoV::Array{Array{Float64}}, X_pölli::Array{Float64}, Y_pölli::Array{Float64}, Z_pölli::Array{Float64}, ρ::Float64, i::Int64, j::Int64)
    if X_pölli[i,j] == X_pölli[i,j+1] ==  X_pölli[i+1,j] == X_pölli[i+1,j+1] && Y_pölli[i,j] == Y_pölli[i,j+1] ==  Y_pölli[i+1,j] == Y_pölli[i+1,j+1]

        I_element = [Ixx 0.0 0.0; 0.0 Iyy 0.0; 0.0 0.0 Izz]

    else
        # Jaetaan elementti 8:aan tetrahedraan, lasketaan tetrahedrojen inertiat elementin keskipisteen suhteen ja summataan yhteen
        Dvec = Array(Float64,3)
        Evec = Array(Float64,3)
        Fvec = Array(Float64,3)
        G = Array(Float64,3)
        H = Array(Float64,3)
        N =  Array(Float64,3)
        volume_tetra::Float64 = 0.0
        Pxx =  Array(Float64,8)
        Pyy =  Array(Float64,8)
        Pzz =  Array(Float64,8)
        Pxy =  Array(Float64,8)
        Pxz =  Array(Float64,8)
        Pyz =  Array(Float64,8)

        origo_temp = CoV[i,j]

        @simd for l in range(1,8)
        # Jaetaan tetrahedroihin
            if l == 1
                @inbounds @fastmath Dvec = [0.0-origo_temp[1]; 0.0-origo_temp[2]; Z_pölli[j]-origo_temp[3]]
                @inbounds @fastmath Evec = [X_pölli[i,j]-origo_temp[1]; Y_pölli[i,j]-origo_temp[2]; Z_pölli[j]-origo_temp[3]]
                @inbounds @fastmath Fvec = [X_pölli[i+1,j]-origo_temp[1]; Y_pölli[i+1,j]-origo_temp[2]; Z_pölli[j]-origo_temp[3]]
            elseif l == 2
                @inbounds @fastmath Dvec = [0.0-origo_temp[1]; 0.0-origo_temp[2]; Z_pölli[j+1]-origo_temp[3]]
                @inbounds @fastmath Evec = [X_pölli[i+1,j+1]-origo_temp[1]; Y_pölli[i+1,j+1]-origo_temp[2]; Z_pölli[j+1]-origo_temp[3]]
                @inbounds @fastmath Fvec = [X_pölli[i,j+1]-origo_temp[1]; Y_pölli[i,j+1]-origo_temp[2]; Z_pölli[j+1]-origo_temp[3]]
            elseif l == 3
                @inbounds @fastmath Dvec = [0.0-origo_temp[1]; 0.0-origo_temp[2]; Z_pölli[j]-origo_temp[3]]
                @inbounds @fastmath Evec = [X_pölli[i+1,j]-origo_temp[1]; Y_pölli[i+1,j]-origo_temp[2]; Z_pölli[j]-origo_temp[3]]
                @inbounds @fastmath Fvec = [0.0-origo_temp[1]; 0.0-origo_temp[2]; Z_pölli[j+1]-origo_temp[3]]
            elseif l == 4
                @inbounds @fastmath Dvec = [X_pölli[i+1,j]-origo_temp[1]; Y_pölli[i+1,j]-origo_temp[2]; Z_pölli[j]-origo_temp[3]]
                @inbounds @fastmath Evec = [X_pölli[i+1,j+1]-origo_temp[1]; Y_pölli[i+1,j+1]-origo_temp[2]; Z_pölli[j+1]-origo_temp[3]]
                @inbounds @fastmath Fvec = [0.0-origo_temp[1]; 0.0-origo_temp[2]; Z_pölli[j+1]-origo_temp[3]]
            elseif l == 5
                @inbounds @fastmath Dvec = [X_pölli[i+1,j]-origo_temp[1]; Y_pölli[i+1,j]-origo_temp[2]; Z_pölli[j]-origo_temp[3]]
                @inbounds @fastmath Evec = [X_pölli[i,j+1]-origo_temp[1]; Y_pölli[i,j+1]-origo_temp[2]; Z_pölli[j+1]-origo_temp[3]]
                @inbounds @fastmath Fvec = [X_pölli[i+1,j+1]-origo_temp[1]; Y_pölli[i+1,j+1]-origo_temp[2]; Z_pölli[j+1]-origo_temp[3]]
            elseif l == 6
                @inbounds @fastmath Dvec = [X_pölli[i+1,j]-origo_temp[1]; Y_pölli[i+1,j]-origo_temp[2]; Z_pölli[j]-origo_temp[3]]
                @inbounds @fastmath Evec = [X_pölli[i,j]-origo_temp[1]; Y_pölli[i,j]-origo_temp[2]; Z_pölli[j]-origo_temp[3]]
                @inbounds @fastmath Fvec = [X_pölli[i,j+1]-origo_temp[1]; Y_pölli[i,j+1]-origo_temp[2]; Z_pölli[j+1]-origo_temp[3]]
            elseif l == 7
                @inbounds @fastmath Dvec = [0.0-origo_temp[1]; 0.0-origo_temp[2]; Z_pölli[j]-origo_temp[3]]
                @inbounds @fastmath Evec = [X_pölli[i,j+1]-origo_temp[1]; Y_pölli[i,j+1]-origo_temp[2]; Z_pölli[j+1]-origo_temp[3]]
                @inbounds @fastmath Fvec = [X_pölli[i,j]-origo_temp[1]; Y_pölli[i,j]-origo_temp[2]; Z_pölli[j]-origo_temp[3]]
            else
                @inbounds @fastmath Dvec = [0.0-origo_temp[1]; 0.0-origo_temp[2]; Z_pölli[j]-origo_temp[3]]
                @inbounds @fastmath Evec = [0.0-origo_temp[1]; 0.0-origo_temp[2]; Z_pölli[j+1]-origo_temp[3]]
                @inbounds @fastmath Fvec = [X_pölli[i,j+1]-origo_temp[1]; Y_pölli[i,j+1]-origo_temp[2]; Z_pölli[j+1]-origo_temp[3]]
            end

            # Lasketaan tetrahedran pohjan reunat
            G = Evec - Dvec
            H = Fvec - Dvec
            # Lasketaan tetrahedran pohjan normaali
            #N[l] = cross(G[l],H[l])
            N = cross(H,G)
            # Lasketaan tetrahedran tilavuus
            volume_tetra = dot(Dvec/3,N/2)

            # Lasketaan tetrahedrojen products of inertiat Pjk, missä j,k = x, y, z = 1, 2, 3
            @inbounds @fastmath Pxx[l] = func_Pjk_tetra2(ρ, volume_tetra, Dvec, Evec, Fvec, 1, 1)
            @inbounds @fastmath Pyy[l] = func_Pjk_tetra2(ρ, volume_tetra, Dvec, Evec, Fvec, 2, 2)
            @inbounds @fastmath Pzz[l] = func_Pjk_tetra2(ρ, volume_tetra, Dvec, Evec, Fvec, 3, 3)

            @inbounds @fastmath Pxy[l] = func_Pjk_tetra2(ρ, volume_tetra, Dvec, Evec, Fvec, 1, 2)
            @inbounds @fastmath Pxz[l] = func_Pjk_tetra2(ρ, volume_tetra, Dvec, Evec, Fvec, 1, 3)
            @inbounds @fastmath Pyz[l] = func_Pjk_tetra2(ρ, volume_tetra, Dvec, Evec, Fvec, 2, 3)
        end

        # Lasketaan elementin inertiat
        Ixx = sum(Pyy) + sum(Pzz)
        Iyy = sum(Pxx) + sum(Pzz)
        Izz = sum(Pxx) + sum(Pyy)
        Iyz = -sum(Pyz)
        Ixz = -sum(Pxz)
        Ixy = -sum(Pxy)

        I_element = [Ixx Ixy Ixz; Ixy Iyy Iyz; Ixz Iyz Izz]
    end

    return I_element
end

function func_inertias_element3(CoV::Array{Float64}, X_pölli11::Float64, X_pölli21::Float64, X_pölli12::Float64, X_pölli22::Float64, Y_pölli11::Float64, Y_pölli21::Float64, Y_pölli12::Float64, Y_pölli22::Float64, Z_pölli1::Float64, Z_pölli2::Float64, ρ::Float64, i::Int64, j::Int64)
    # Jaetaan elementti 8:aan tetrahedraan, lasketaan tetrahedrojen inertiat elementin keskipisteen suhteen ja summataan yhteen
    Dvec = Array(Float64,3)
    Evec = Array(Float64,3)
    Fvec = Array(Float64,3)
    G = Array(Float64,3)
    H = Array(Float64,3)
    N =  Array(Float64,3)
    volume_tetra::Float64 = 0.0
    Pxx =  Array(Float64,8)
    Pyy =  Array(Float64,8)
    Pzz =  Array(Float64,8)
    Pxy =  Array(Float64,8)
    Pxz =  Array(Float64,8)
    Pyz =  Array(Float64,8)

    #origo_temp = CoV

    @simd for l in range(1,8)
    # Jaetaan tetrahedroihin
        if l == 1
            @inbounds @fastmath Dvec = [0.0-CoV[1]; 0.0-CoV[2]; Z_pölli1-CoV[3]]
            @inbounds @fastmath Evec = [X_pölli11-CoV[1]; Y_pölli11-CoV[2]; Z_pölli1-CoV[3]]
            @inbounds @fastmath Fvec = [X_pölli21-CoV[1]; Y_pölli21-CoV[2]; Z_pölli1-CoV[3]]
        elseif l == 2
            @inbounds @fastmath Dvec = [0.0-CoV[1]; 0.0-CoV[2]; Z_pölli2-CoV[3]]
            @inbounds @fastmath Evec = [X_pölli22-CoV[1]; Y_pölli22-CoV[2]; Z_pölli2-CoV[3]]
            @inbounds @fastmath Fvec = [X_pölli12-CoV[1]; Y_pölli12-CoV[2]; Z_pölli2-CoV[3]]
        elseif l == 3
            @inbounds @fastmath Dvec = [0.0-CoV[1]; 0.0-CoV[2]; Z_pölli1-CoV[3]]
            @inbounds @fastmath Evec = [X_pölli21-CoV[1]; Y_pölli21-CoV[2]; Z_pölli1-CoV[3]]
            @inbounds @fastmath Fvec = [0.0-CoV[1]; 0.0-CoV[2]; Z_pölli2-CoV[3]]
        elseif l == 4
            @inbounds @fastmath Dvec = [X_pölli21-CoV[1]; Y_pölli21-CoV[2]; Z_pölli1-CoV[3]]
            @inbounds @fastmath Evec = [X_pölli22-CoV[1]; Y_pölli22-CoV[2]; Z_pölli2-CoV[3]]
            @inbounds @fastmath Fvec = [0.0-CoV[1]; 0.0-CoV[2]; Z_pölli2-CoV[3]]
        elseif l == 5
            @inbounds @fastmath Dvec = [X_pölli21-CoV[1]; Y_pölli21-CoV[2]; Z_pölli1-CoV[3]]
            @inbounds @fastmath Evec = [X_pölli12-CoV[1]; Y_pölli12-CoV[2]; Z_pölli2-CoV[3]]
            @inbounds @fastmath Fvec = [X_pölli22-CoV[1]; Y_pölli22-CoV[2]; Z_pölli2-CoV[3]]
        elseif l == 6
            @inbounds @fastmath Dvec = [X_pölli21-CoV[1]; Y_pölli21-CoV[2]; Z_pölli1-CoV[3]]
            @inbounds @fastmath Evec = [X_pölli11-CoV[1]; Y_pölli11-CoV[2]; Z_pölli1-CoV[3]]
            @inbounds @fastmath Fvec = [X_pölli12-CoV[1]; Y_pölli12-CoV[2]; Z_pölli2-CoV[3]]
        elseif l == 7
            @inbounds @fastmath Dvec = [0.0-CoV[1]; 0.0-CoV[2]; Z_pölli1-CoV[3]]
            @inbounds @fastmath Evec = [X_pölli12-CoV[1]; Y_pölli12-CoV[2]; Z_pölli2-CoV[3]]
            @inbounds @fastmath Fvec = [X_pölli11-CoV[1]; Y_pölli11-CoV[2]; Z_pölli1-CoV[3]]
        else
            @inbounds @fastmath Dvec = [0.0-CoV[1]; 0.0-CoV[2]; Z_pölli1-CoV[3]]
            @inbounds @fastmath Evec = [0.0-CoV[1]; 0.0-CoV[2]; Z_pölli2-CoV[3]]
            @inbounds @fastmath Fvec = [X_pölli12-CoV[1]; Y_pölli12-CoV[2]; Z_pölli2-CoV[3]]
        end

        # Lasketaan tetrahedran pohjan reunat
        G = Evec - Dvec
        H = Fvec - Dvec
        # Lasketaan tetrahedran pohjan normaali
        #N[l] = cross(G[l],H[l])
        N = cross(H,G)
        # Lasketaan tetrahedran tilavuus
        volume_tetra = dot(Dvec/3,N/2)

        # Lasketaan tetrahedrojen products of inertiat Pjk, missä j,k = x, y, z = 1, 2, 3
        @inbounds @fastmath Pxx[l] = func_Pjk_tetra2(ρ, volume_tetra, Dvec, Evec, Fvec, 1, 1)
        @inbounds @fastmath Pyy[l] = func_Pjk_tetra2(ρ, volume_tetra, Dvec, Evec, Fvec, 2, 2)
        @inbounds @fastmath Pzz[l] = func_Pjk_tetra2(ρ, volume_tetra, Dvec, Evec, Fvec, 3, 3)

        @inbounds @fastmath Pxy[l] = func_Pjk_tetra2(ρ, volume_tetra, Dvec, Evec, Fvec, 1, 2)
        @inbounds @fastmath Pxz[l] = func_Pjk_tetra2(ρ, volume_tetra, Dvec, Evec, Fvec, 1, 3)
        @inbounds @fastmath Pyz[l] = func_Pjk_tetra2(ρ, volume_tetra, Dvec, Evec, Fvec, 2, 3)
    end

    # Lasketaan elementin inertiat
    Ixx = sum(Pyy) + sum(Pzz)
    Iyy = sum(Pxx) + sum(Pzz)
    Izz = sum(Pxx) + sum(Pyy)
    Iyz = -sum(Pyz)
    Ixz = -sum(Pxz)
    Ixy = -sum(Pxy)

    I_element = [Ixx Ixy Ixz; Ixy Iyy Iyz; Ixz Iyz Izz]

    return I_element
end



# Laskee tetrahedran product of inertian Pjk:n, missä j,k = x, y, z
function func_Pjk_tetra(ρ::Float64, volume_tetra::Array{Float64}, Dvec::Array{Array{Float64}}, Evec::Array{Array{Float64}}, Fvec::Array{Array{Float64}}, i::Int64, j::Int64, k::Int64)
    Pjk_tetra = (ρ*volume_tetra[i]/20.0) * (2*Dvec[i][j]*Dvec[i][k] + 2*Evec[i][j]*Evec[i][k] + 2*Fvec[i][j]*Fvec[i][k] + Dvec[i][j]*Evec[i][k] + Dvec[i][k]*Evec[i][j] + Dvec[i][j]*Fvec[i][k] + Dvec[i][k]*Fvec[i][j] + Evec[i][j]*Fvec[i][k] + Evec[i][k]*Fvec[i][j])
end

function func_Pjk_tetra2(ρ::Float64, volume_tetra::Float64, Dvec::Array{Float64}, Evec::Array{Float64}, Fvec::Array{Float64}, j::Int64, k::Int64)
    Pjk_tetra = (ρ*volume_tetra/20.0) * (2*Dvec[j]*Dvec[k] + 2*Evec[j]*Evec[k] + 2*Fvec[j]*Fvec[k] + Dvec[j]*Evec[k] + Dvec[k]*Evec[j] + Dvec[j]*Fvec[k] + Dvec[k]*Fvec[j] + Evec[j]*Fvec[k] + Evec[k]*Fvec[j])
end


# Laskee pöllin kaikkien elementtien inertiat ja tallentaa ne matriisiin.
function func_inertias_elements(CoV::Array{Array{Float64}}, X_pölli::Array{Float64}, Y_pölli::Array{Float64}, Z_pölli::Array{Float64}, ρ::Float64, n::Int64, k::Int64)
    I_elements = Array(Array{Float64},n,k-1)

    # Lasketaan elementtien inertiat
    @simd for j in range(1,k-1)
        @simd for i in range(1,n)
#            @inbounds @fastmath I_elements[i,j] = func_inertias_element(CoV, X_pölli, Y_pölli, Z_pölli, ρ, i, j)
            @inbounds @fastmath I_elements[i,j] = func_inertias_element2(CoV, X_pölli, Y_pölli, Z_pölli, ρ, i, j)
#            @inbounds @fastmath I_elements[i,j] = func_inertias_element3(CoV[i,j], X_pölli[i,j], X_pölli[i+1,j], X_pölli[i,j+1], X_pölli[i+1,j+1], Y_pölli[i,j], Y_pölli[i+1,j], Y_pölli[i,j+1], Y_pölli[i+1,j+1], Z_pölli[j], Z_pölli[j+1], ρ::Float64, i::Int64, j::Int64)
        end
    end

return I_elements
end

# Laskee pöllin inertian pöllin keskipisteen suhteen
function func_inertias_pölli(CoV::Array{Array{Float64}}, I_elements::Array{Array{Float64}}, volume_elements::Array{Float64}, ρ::Float64, x_CoV::Float64, y_CoV::Float64, z_CoV::Float64, n::Int64, k::Int64)
    # Muutetaan elementtien inertiat pöllin keskipisteen suhteen
    # Parallel axis theorem
    I_elements_translation = deepcopy(I_elements)
    I_total = zeros(3,3)

    @simd for j in range(1,k-1)
        @simd for i in range(1,n)
            # First column
            # Ixx
            @inbounds @fastmath I_elements_translation[i,j][1,1] = I_elements_translation[i,j][1,1] + ρ*volume_elements[i,j]*( (CoV[i,j][2]-y_CoV)^2 + (CoV[i,j][3]-z_CoV)^2 )
            # Ixy
            @inbounds @fastmath I_elements_translation[i,j][2,1] = I_elements_translation[i,j][2,1] + ρ*volume_elements[i,j]*( (CoV[i,j][1]-x_CoV) * (CoV[i,j][2]-y_CoV))
            # Ixz
            @inbounds @fastmath I_elements_translation[i,j][3,1] = I_elements_translation[i,j][3,1] + ρ*volume_elements[i,j]*( (CoV[i,j][1]-x_CoV) * (CoV[i,j][3]-z_CoV))
            # Second column
            # Ixy
            @inbounds @fastmath I_elements_translation[i,j][1,2] = I_elements_translation[i,j][2,1]
            # Iyy
            @inbounds @fastmath I_elements_translation[i,j][2,2] = I_elements_translation[i,j][2,2] + ρ*volume_elements[i,j]*( (CoV[i,j][1]-x_CoV)^2 + (CoV[i,j][3]-z_CoV)^2 )
            # Iyz
            @inbounds @fastmath I_elements_translation[i,j][3,2] = I_elements_translation[i,j][3,2] + ρ*volume_elements[i,j]*( (CoV[i,j][2]-y_CoV) * (CoV[i,j][3]-z_CoV))
            # Third column
            # Ixz
            @inbounds @fastmath I_elements_translation[i,j][1,3] = I_elements_translation[i,j][3,1]
            # Iyz
            @inbounds @fastmath I_elements_translation[i,j][2,3] = I_elements_translation[i,j][3,2]
            # Izz
            @inbounds @fastmath I_elements_translation[i,j][3,3] = I_elements_translation[i,j][3,3] + ρ*volume_elements[i,j]*( (CoV[i,j][1]-x_CoV)^2 + (CoV[i,j][2]-y_CoV)^2 )
        end
    end

    # Summataan elementtien inertiakomponentit keskenään.
    @simd for j in range(1,k-1)
        @simd for i in range(1,n)
            @inbounds @fastmath I_total[1,1] = I_total[1,1] + I_elements_translation[i,j][1,1]
            @inbounds @fastmath I_total[2,1] = I_total[2,1] + I_elements_translation[i,j][2,1]
            @inbounds @fastmath I_total[3,1] = I_total[3,1] + I_elements_translation[i,j][3,1]

            @inbounds @fastmath I_total[2,2] = I_total[2,2] + I_elements_translation[i,j][2,2]
            @inbounds @fastmath I_total[3,2] = I_total[3,2] + I_elements_translation[i,j][3,2]

            @inbounds @fastmath I_total[3,3] = I_total[3,3] + I_elements_translation[i,j][3,3]

        end
    end
    I_total[1,2] = I_total[2,1]
    I_total[1,3] = I_total[3,1]
    I_total[2,3] = I_total[3,2]

return I_total
end

#Testausta varten
R_sektori_1 = [0.116; 0.121; 0.114; 0.119; 0.122; 0.122; 0.115; 0.107; 0.126; 0.116]
R_sektori_2 = [0.116; 0.113; 0.126; 0.12; 0.125; 0.121; 0.126; 0.115; 0.115; 0.116]
R_sektori_3 = [0.119; 0.109; 0.121; 0.11; 0.119; 0.121; 0.115; 0.115; 0.109; 0.119]

R_sektori = Array(Float64,10,3)
R_sektori[:,1] = R_sektori_1
R_sektori[:,2] = R_sektori_2
R_sektori[:,3] = R_sektori_3

Θ_sektori = [
 0.0
 0.698132
 1.39626
 2.0944
 2.79253
 3.49066
 4.18879
 4.88692
 5.58505
 6.28319]