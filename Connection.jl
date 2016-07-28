
println("Waiting for connection")
server = listen(5111) # Avataan Modelerissa asetettu portti

conn = accept(server)

ninputs = 23 # Mevea solverista tulevat outputit
noutputs = 0  # Mevea solveriin menevät inputit

params=Array{Int32}(3)
params2=Array{Int32}(3)
params[1]=1
params[2]=noutputs
params[3]=ninputs
params2 = read(conn,Int32,3) #Tarksitetaan että Ratkaisiassa ja Juliassa on sama määrä out-ja inputteja.

write(conn,params)

ins=Array{Float64}(ninputs)
outs=Array{Float64}(noutputs)
All=Array{Float64}(ninputs)
Rotation_matrix_1=Array{Float64}(3,3)
time_step=0.001
time=0

println("Connection established")
while isopen(conn) #Kommunikointi

  try
    ins = read(conn,Float64,ninputs)
    append!(All,ins)
    #ins[1] Simulation time [s]
    #ins[2] Knife bar blade angle [deg(90,-90)]
    #ins[3] Knife bar blade center Global location x [m]
    #ins[4] Knife bar blade center Global location y [m]
    #ins[5] Log global location x [m]
    #ins[6] Log global location y [m]
    #ins[7] Log global location z [m]
    #ins[8] Log euler parameter e0
    #ins[9] Log euler parameter e1
    #ins[10] Log euler parameter e2
    #ins[11] Log euler parameter e3
    #ins[12] Lathe spindle R Global location z [m]
    #ins[13] Lathe spindle L Global location z [m]
    #ins[14] Lathe spindle R orientation z [deg[90,-90]]
    #ins[15] Lathe spindle L orientation z [deg[90,-90]]
    #ins[16] Lathe spindle R angular velocity z [rad/s]
    #ins[17] Lathe spindle L angular velocity z [rad/s]
    #ins[18] Round bar center Global location x [m]
    #ins[19] Round bar center Global location y [m]
    #ins[20] Back up roll1 Global location x [m]
    #ins[21] Back up roll1 Global location y [m]
    #ins[22] Back up roll2 Global location x [m]
    #ins[23] Back up roll2 Global location y [m]
    e0=ins[8]
    e1=ins[9]
    e2=ins[10]
    e3=ins[11]
    
    Rotation_matrix_1[1,1]=1 - 2 * e2^2 - 2 * e3^2
    Rotation_matrix_1[1,2]=2*(e1 * e2 - e0 * e3)
    Rotation_matrix_1[1,3]=2*(e1 * e3 + e0 * e2)
    
    Rotation_matrix_1[2,1]=2*(e1 * e2 + e0 * e3)
    Rotation_matrix_1[2,2]=1 - 2 * e1^2 - 2 * e3^2
    Rotation_matrix_1[2,3]=2*(e2 * e3 - e0 * e1)
    
    Rotation_matrix_1[3,1]=2*(e1 * e3 - e0 * e2)
    Rotation_matrix_1[3,2]=2*(e2 * e3 + e0 * e1)
    Rotation_matrix_1[3,3]=1 - 2 * e1^2 - 2 * e2^2
    write(conn,outs)
  catch
    break
  end
end

println("Connection closed")
close(conn)
close(server)
