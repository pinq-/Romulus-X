# Lisätään Julialle koneen ytimet käyttöön
# addprocs(CPU_CORES - 1)

####################
# TIME IT!
function timeit()
    @everywhere k = 486
    @everywhere n = 360
    @everywhere R_alku = 0.12
    a = Array(Float64,n,k)
    q = SharedArray(Float64, n, k)

    @time advection_shared!(q)
    @time testi_p(k::Int64,n::Int64,a)
    #println(q == a)
end


# This function retuns the (irange,jrange) indexes assigned to this worker
# It splits the q::SharedArray into chunks of columns for workers
@everywhere function myrange(q::SharedArray)
    idx = indexpids(q)
    if idx == 0
        # This worker is not assigned a piece
        return 1:0, 1:0
    end
    nchunks = length(procs(q))
    splits = [round(Int, s) for s in linspace(0,size(q,2),nchunks+1)]
    1:size(q,1), splits[idx]+1:splits[idx+1]
end

# This function retuns the (irange,jrange) indexes assigned to this worker
# It splits the q::SharedArray into chunks of rows for workers
@everywhere function myrangerows(q::SharedArray)
    idx = indexpids(q)
    if idx == 0
        # This worker is not assigned a piece
        return 1:0, 1:0
    end
    nchunks = length(procs(q))
    splits = [round(Int, s) for s in linspace(0,size(q,1),nchunks+1)]
    splits[idx]+1:splits[idx+1], 1:size(q,2)
end


# Here's the kernel
@everywhere function advection_chunk!(q, irange, jrange)
    @show (irange, jrange)  # display so we can see what's happening
    for j in jrange, i in irange
        q[i,j] = R_alku * (rand(89:105)/100)
    end
    q
end

# Here's a convenience wrapper for a SharedArray implementation
@everywhere advection_shared_chunk!(q) = advection_chunk!(q, myrange(q)...)

function advection_shared!(q)
    @sync begin
        for p in procs(q)
            @async remotecall_wait(advection_shared_chunk!, p, q)
        end
    end
    q
end

@everywhere function testi_p(k::Int64,n::Int64,a_shared)
    @simd for j in range(1,k)
        @simd for i in range(1,n)
            @fastmath @inbounds a_shared[i,j] = R_alku * (rand(89:105)/100)
        end
    end
    return a_shared
end
