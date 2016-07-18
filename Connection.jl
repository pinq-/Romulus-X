println("Waiting for connection")
server = listen(5111) # Avataan Modelerissa asetettu portti

conn = accept(server)

ninputs = 3 # Mevea solverista tulevien outputien määrä
noutputs = 2  # Mevea solveriin menevien inputien määrä

params=Array{Int32}(3) 
params2=Array{Int32}(3)
params[1]=1
params[2]=noutputs 
params[3]=ninputs
params2 = read(conn,Int32,3) #Tarksitetaan että ratkaisiassa ja Juliassa on sama määrä out-ja inputteja.
write(conn,params)

ins=Array{Float64}(ninputs)
outs=Array{Float64}(noutputs)

println(")Connection established")
while isopen(conn) #Kommunikointi

  try
    ins = read(conn,Float64,ninputs)
    println(ins)
    write(conn,outs)
  catch
    break
  end
end
println("Connection closed")
close(conn)
close(server)

