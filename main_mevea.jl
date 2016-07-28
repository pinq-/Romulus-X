include("lahto_arvot.jl")
include("3D_laskentafunktiot.jl")
include("apufunktiot.jl")

using PyCall
using PyPlot
@pyimport matplotlib.cm as cm


function connection_Mevea()
    #################################
    # Alustus
    k::Int64
    k = 3
    n::Int64
    n = 179
    w::Float64
    w = 1.35

    ##############
    # Terä
    X_terä = Array(Float64, k)
    Y_terä = Array(Float64, k)
    Z_terä = Array(Float64, k)

    # Lasketaan terän pituus
    Z_terä[1] = -w/2
    Z_terä[k] = w/2
    w_terä = w / (length(Z_terä)-1)
    @simd for i in range(2,k-1)
        #if i == 1
            #Z_terä[i] = ins[7]
        #elseif i < k
            @fastmath @inbounds Z_terä[i] = Z_terä[1] + w_terä*(i-1)
            #@fastmath @inbounds tera_asema[i] = tera_asema[1] + (Z_terä[i] * tan(tera_apukulma))
        #else
            #@fastmath @inbounds Z_terä[i] = ins[10]
        #end
    end
    ##############
    # Pöllin lokaalit koordinaatit
    X_pölli = Array(Float64, n+1, k)
    Y_pölli = Array(Float64, n+1, k)
    Z_pölli = func_z_serial_mevea(k,w)
    XYZ_pölli_local = Array(Tuple{Float64,Float64,Float64}, n+1, k)
    # Pöllin globaalit koordinaatit
    #X_pölli_global = Array(Float64, n+1, k)
    #Y_pölli_global = Array(Float64, n+1, k)
    #Z_pölli_global = Array(Float64, k)
    XYZ_pölli_global = Array(Tuple{Float64,Float64,Float64}, n+1, k)
    XYZ_pölli_global_temp = Array(Float64, 3)
    # Muodostetaan sylinterimäinen pölli polaarikoordinaattien avulla
    Θ_sektori, R_sektori = func_sektorit_serial(n,k,R_alku)
    # Muutetaan polaarikoordinaatit takaisin pöllin lokaaliin karteesiseen koordinaatistoon
    func_XY_pölli_serial!(X_pölli,Y_pölli,R_sektori, Θ_sektori, n, k)
    # Luodaan n,k kokoinen Array of tuples
    for j in range(1,k)
        for i in range(1,n+1)
            XYZ_pölli_local[i,j] = (X_pölli[i,j],Y_pölli[i,j],Z_pölli[j])
        end
    end

    println(XYZ_pölli_local[1,1])

    #################################
    # Määrittää kuvaajan - rotaatiomatriisin elementeille
    fig, axs = plt[:subplots](3, 3)

    for i in range(1,3)
        for j in range(1,3)
            axs[i,j][:set_title]("[$(i),$(j)]")
            axs[i,j][:set_xlim](0.0,300.0)
            axs[i,j][:set_ylim](-2.0,2.2)
        end
    end

    R_matrix = Array(Array{Float64},3,3)
    tempvector = Array(Float64,0)
    for i in range(1,3)
        for j in range(1,3)
            R_matrix[i,j] = copy(tempvector)
        end
    end

    aika = Array(Float64,0)

    #################################
    # Yhteys Meveaan
    #################################
    println("Waiting for connection")
    server = listen(5111) # Avataan Modelerissa asetettu portti

    conn = accept(server)

    ninputs = 29 # Mevea solverista tulevat outputit
    noutputs = 0  # Mevea solveriin menevät inputit

    params=Array{Int32}(3)
    params2=Array{Int32}(3)
    params[1]=1
    params[2]=noutputs
    params[3]=ninputs
    params2 = read(conn,Int32,3) #Tarkistetaan että Ratkaisiassa ja Juliassa on sama määrä out-ja inputteja.

    write(conn,params)

    ins=Array{Float64}(ninputs)
    outs=Array{Float64}(noutputs)
    All=Array{Float64}(ninputs)
    Rotation_matrix_1=Array{Float64}(3,3)
    time_step=0.001
    time=0

#=    Rotation_matrix_1[1,1]=ins[14]^2 + ins[15]^2 - ins[16]^2 - ins[17]^2
    Rotation_matrix_1[1,2]=2*(ins[15]*ins[16] - ins[14]*ins[17])
    Rotation_matrix_1[1,3]=2*(ins[15]*ins[17] + ins[14]*ins[16])
    Rotation_matrix_1[2,1]=2*(ins[15]*ins[16] + ins[14]*ins[17])
    Rotation_matrix_1[2,2]=ins[14]^2 - ins[15]^2 + ins[16]^2 - ins[17]^2
    Rotation_matrix_1[2,3]=2*(ins[16]*ins[17] - ins[14]*ins[15])
    Rotation_matrix_1[3,1]=2*(ins[15]*ins[17] - ins[14]*ins[16])
    Rotation_matrix_1[3,2]=2*(ins[16]*ins[17] + ins[14]*ins[15])
    Rotation_matrix_1[3,3]=ins[14]^2 - ins[15]^2 - ins[16]^2 + ins[17]^2
=#

    ax1, = axs[1,1][:plot](aika,R_matrix[1,1])
    ax2, = axs[1,2][:plot](aika,R_matrix[1,2])
    ax3, = axs[1,3][:plot](aika,R_matrix[1,3])
    ax4, = axs[2,1][:plot](aika,R_matrix[2,1])
    ax5, = axs[2,2][:plot](aika,R_matrix[2,2])
    ax6, = axs[2,3][:plot](aika,R_matrix[2,3])
    ax7, = axs[3,1][:plot](aika,R_matrix[3,1])
    ax8, = axs[3,2][:plot](aika,R_matrix[3,2])
    ax9, = axs[3,3][:plot](aika,R_matrix[3,3])


    println("Connection established")
    while isopen(conn) #Kommunikointi

      try
        ins = read(conn,Float64,ninputs)
        #append!(All,ins)
        #ins[1] Simulation time [s]
        #ins[2] Knife bar blade angle [deg(90,-90)]
        #ins[3] Knife bar blade center Global location x [m]
        #ins[4] Knife bar blade center Global location y [m]
        #ins[5] Knife bar blade left Global location x [m]
        #ins[6] Knife bar blade left Global location y [m]
        #ins[7] Knife bar blade left Global location z [m]
        #ins[8] Knife bar blade right Global location x [m]
        #ins[9] Knife bar blade right Global location y [m]
        #ins[10] Knife bar blade right Global location z [m]
        #ins[11] Log global location x [m]
        #ins[12] Log global location y [m]
        #ins[13] Log global location z [m]
        #ins[14] Log euler parameter e0
        #ins[15] Log euler parameter e1
        #ins[16] Log euler parameter e2
        #ins[17] Log euler parameter e3
        #ins[18] Lathe spindle R Global location z [m]
        #ins[19] Lathe spindle L Global location z [m]
        #ins[20] Lathe spindle R orientation z [deg[90,-90]]
        #ins[21] Lathe spindle L orientation z [deg[90,-90]]
        #ins[22] Lathe spindle R angular velocity z [rad/s]
        #ins[23] Lathe spindle L angular velocity z [rad/s]
        #ins[24] Round bar center Global location x [m]
        #ins[25] Round bar center Global location y [m]
        #ins[26] Back up roll1 Global location x [m]
        #ins[27] Back up roll1 Global location y [m]
        #ins[28] Back up roll2 Global location x [m]
        #ins[29] Back up roll2 Global location y [m]

        # Pöllin rotaatiomatriisi
        #println("Timing")
        Rotation_matrix_1 = func_rotation_matrix(ins[14],ins[15],ins[16],ins[17])
        for j in range(1,3)
            for i in range(1,3)
                push!(R_matrix[i,j],Rotation_matrix_1[i,j])
            end
        end
        push!(aika,ins[1])


        # Pöllin keskipisteen asema globaalissa koordinaatistossa
        XYZ_log = (ins[11],ins[12],ins[13])
        # Pöllin profiilien koordinaatit globaalissa koordinaatistossa
        for j in range(1,k)
            @simd for i in range(1,n+1)
                @inbounds @fastmath XYZ_pölli_global_temp = [XYZ_log[1];XYZ_log[2];XYZ_log[3]] + Rotation_matrix_1*[XYZ_pölli_local[i,j][1];XYZ_pölli_local[i,j][2];XYZ_pölli_local[i,j][3]]
                @inbounds XYZ_pölli_global[i,j] = (XYZ_pölli_global_temp[1],XYZ_pölli_global_temp[2],XYZ_pölli_global_temp[3])
            end
        end
        #X_pölli_global = XYZ_log + Rotation_matrix_1*X_pölli
        #Y_pölli_global = XYZ_log + Rotation_matrix_1*Y_pölli
        #Z_pölli_global = XYZ_log + Rotation_matrix_1*Z_pölli

#=        X_terä[1] = ins[5]
        Y_terä[1] = ins[6]
        X_terä[k] = ins[8]
        Y_terä[k] = ins[9]=#


        #println("XYZ_log = ", XYZ_log)
        #println(XYZ_pölli_global[1,1])
        #println(XYZ_pölli_global[1,2])
        #println(XYZ_pölli_global[1,3])

        #p = [ins[14]; ins[15]; ins[16]; ins[17]]

        #testi = p.' * p - 1
        #if testi != 0.0
        #    println("Eulerin parametrit vituillaan!")
        #    println(testi)
        #end
        #println(Rotation_matrix_1)
        #println("e1 = ", ins[14])
        #println("e2 = ", ins[15])
        #println("e3 = ", ins[16])
        #println("e4 = ", ins[17])









        write(conn,outs)
      catch
        break
      end
    end

    println("Connection closed")
    close(conn)
    close(server)

    ax1[:set_data](aika,R_matrix[1,1])
    ax2[:set_data](aika,R_matrix[1,2])
    ax3[:set_data](aika,R_matrix[1,3])
    ax4[:set_data](aika,R_matrix[2,1])
    ax5[:set_data](aika,R_matrix[2,2])
    ax6[:set_data](aika,R_matrix[2,3])
    ax7[:set_data](aika,R_matrix[3,1])
    ax8[:set_data](aika,R_matrix[3,2])
    ax9[:set_data](aika,R_matrix[3,3])

end #function end
