# makes calling C functions a bit easier
macro ccd_ccall(func, args...)
    appfolder = dirname(@__FILE__)
    lib = "$appfolder/libconcorde.so" # path to concorde shared lib
    args = map(esc,args)
    @static if Sys.isunix()
        return quote
            ccall(($func,$lib), $(args...))
        end
    end
    @static if Sys.iswindows()
        return quote
            ccall(($func,$lib), $(esc(:stdcall)), $(args...))
        end
    end
end

# function for computing euclydian distance
#=
function ccd_distance(data::DataCluCVRP, arc::Tuple{Int64, Int64})
    e = (arc[1] < arc[2]) ? arc : (arc[2],arc[1])
    u, v = arc
    vertices = data.G′.V′ 
    # array <vertices> is indexed from 1 (depot is vertices[1], customer 1 is vertices[2], and so on) 
    x_sq = (vertices[v].pos_x - vertices[u].pos_x)^2
    y_sq = (vertices[v].pos_y - vertices[u].pos_y)^2
    if data.round
        return floor(sqrt(x_sq + y_sq) + 0.5)
    end
    return sqrt(x_sq + y_sq)
end
=#

# struct required by concorde itself
mutable struct rand_state
    a::Cint
    b::Cint
    arr::NTuple{55, Cint}    
    rand_state() = new(0,0,ntuple(x->0,55))
end

function concorde_st(n,cost)
    status = 0 # Concorde functions return 1 if something fails
    success = Ref(Int32(-1)) # 1 if the run finished normally, and set to 0 if the search was terminated early (by hitting some predefined limit)
    found_tour = Ref(Int32(-1)) # 1 if a tour has been found (if success is 0, then it may not be the optimal tour)
    solving_time = Float64(0.0)
    time_bound = C_NULL # Run time limit (it can be NULL)
    hit_timebound = Ref(Int32(0)) # 1 if timebound was reached
    opt_val = Ref(Float64(0.0)) 
    in_val = C_NULL # Can be used to specify an initial upperbound (it can be NULL)
    
    print_level = 1 # Suppress most output if set to a nonzero value

    # Initialize the portable random number generator
    seed = rand(Int32)
    rstate = rand_state() 
    @ccd_ccall(:CCutil_sprand, Cvoid, (Cint,Ref{rand_state}), seed, rstate) 

    ncount = n # Number of nodes
    ecount = (ncount * (ncount - 1)) / 2 # Number of edges
    
    #nodes = cluster_customers(data,k) # Original nodes
    M = 1e4
    edge = 1
    edgeWeight = 1
    elist = Vector{Cint}(undef, Int64(ecount*2)) # Array giving the ends of the edges (in pairs)
    elen = Vector{Cint}(undef, Int64(ecount)) # Array giving the weights of the edges
    for i in 0:ncount-1 # Concorde considers nodes ids starting at 0
        for j in i+1:ncount-1
            if i != j
                elist[edge] = i;
                elist[edge + 1] = j;

                #elen[edgeWeight] = ccd_distance(data, (nodes[i+1],nodes[j+1])) + M
                elen[edgeWeight] = cost[i+1,j+1] + M
                
                edgeWeight += 1;
                edge += 2;
            end
        end
    end

    # Specifes a char string that will be used to name various files that are written during the branch and bound search
    name = @ccd_ccall(:CCtsp_problabel, Cstring, (Cstring,), "_")
    # @show unsafe_string(name)
    
    # Gives a starting tour in node node node format (it can be NULL)
    in_tour = Ptr{Cint}(C_NULL) 
    # Optimal tour (it can be NULL, if it is not NULL then it should point to an array of length at least ncount)
    out_tour = Vector{Cint}(undef, Int64(ncount)) 

    # Start time
    solving_time = @ccd_ccall(:CCutil_zeit, Cdouble, ())

    # Solve TSP
    appfolder = dirname(@__FILE__)
    open("$(appfolder)/out/concorde.log", "w") do fileout
        redirect_stdout(fileout) do
            redirect_stderr(fileout) do
                flush(stdout)            
                status = @ccd_ccall(:CCtsp_solve_sparse, Cint, 
                (Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ref{Cint}, Ptr{Cdouble}, Ref{Cdouble}, 
                Ref{Cint}, Ref{Cint}, Cstring, Ptr{Cdouble}, Ref{Cint}, Cint, Ref{rand_state}),    
                ncount, ecount, elist, elen, in_tour, out_tour, in_val, opt_val, 
                success, found_tour, name, time_bound, hit_timebound, print_level, rstate)
                Base.Libc.flush_cstdio()
                flush(stdout)
            end
        end
    end
    # End time
    solving_time = @ccd_ccall(:CCutil_zeit, Cdouble, ()) - solving_time
    opt_val[] -= (ncount - 1) * M # Solution cost

    (status == 1) && error("Something failed on Concorde functions!")
    (found_tour[] == 0) && error("No tour has been found!")
    if (success[] == 0) && (found_tour[] == 1)
        error("The tour may not be the optimal tour!")
    end

    # Get path with original node ids
    aux = []
    pos_s, pos_t = 0, 0

    tour = []
    for (i,v) in enumerate(out_tour)
         push!(tour,v)
    end

    cost_ = 0.0
    for i=1:length(tour)-1
        cost_ += cost[tour[i]+1,tour[i+1]+1]
    end
    cost_ += cost[tour[end]+1,tour[1]+1]


    #=
    path = []
    cost = 0
    if !(get_path)
        path = aux
        cost = ccd_distance(data, (path[1],path[end]))
        cost -= ccd_distance(data, (path[pos_s],path[pos_t]))
    else
        # It reorders path beginning at s and finishing at t
        if pos_s == 1 && pos_t == ncount
            path = aux
        elseif pos_s == ncount && pos_t == 1
            path = reverse(aux)
        else
            path = (pos_s < pos_t) ? reverse(vcat(aux[pos_t:end],aux[1:pos_s])) : vcat(aux[pos_s:end],aux[1:pos_t])
        end
    end
    # It checks solution cost
    for i in 1:ncount-1
        cost += ccd_distance(data, (path[i],path[i+1]))
    end
    (opt_val[] != cost) && error("The calculated path cost ($cost) does not match with the Concorde solution cost ($(opt_val[]))!")
    
    # @show solving_time, Int64(cost), Int64(opt_val[]), path
    return cost, path
    =#
    return cost_, tour
end

#=
function enumerate_tsp(data::DataCluCVRP, k::Int64, s::Int64, t::Int64)
    if cluster_size(data,k) == 4
        r = setdiff(cluster_customers(data,k),[s,t])
        cost = ccd_distance(data,(s,r[1])) + ccd_distance(data,(r[2],t))
        if cost > ccd_distance(data,(s,r[2])) + ccd_distance(data,(r[1],t))
            cost = ccd_distance(data,(s,r[2])) + ccd_distance(data,(r[1],t))
            r[1], r[2] = r[2], r[1]
        end
        cost += ccd_distance(data,(r[1],r[2]))
        path = [s,r[1],r[2],t]
    elseif cluster_size(data,k) == 3
        r = setdiff(cluster_customers(data,k),[s,t])
        cost = ccd_distance(data,(s,r[1])) + ccd_distance(data,(r[1],t))
        path = [s,r[1],t]  
    else
        cost = ccd_distance(data,(s,t))
        path = [s,t]    
    end
    return cost, path
end

function concorde(data::DataCluCVRP, k::Int64, s::Int64, t::Int64; get_path::Bool=false)
    if cluster_size(data,k) > 4 # concorde olny works for n > 4
    	# This version considers that a path must start at node 's' and finishes at node 't'
    	# To eliinate this condition small changes are required in concorde_st(data, k, s, t, get_path) 
        cost, path = concorde_st(data, k, s, t, get_path)
    else # for n <= 4 all solututions are enumerated
        cost, path = enumerate_tsp(data, k, s, t)
    end
    
    if !(get_path)
        return cost
    else
        return path
    end
end
=#
