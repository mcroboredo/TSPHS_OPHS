#= Remeber that the DataOPHS struct has the following fields:
    - G′::InputGraph
    - Hotels::Array{Int64} # Hotel nodes. The first two hotels are respectively the start depot and the end depot 
    - Customers::Array{Int64} # Customer nodes
    - Tmax::Float64 # Total tour length
    - D::Int64 # Number of trips
    - Td::Array{Float64} # array containing the length of each trip
    - symmetric::Bool # = true if Td[1] = Td[2] = ... = Td[D]
    - cuts::Int64 #Total number of cuts inserted by the user callback
=#


# Está sendo utilizado? Parece que não
mutable struct Layer
    type::String
    l_number::Int64 
    first_id::Int64  
    last_id::Int64 
end


# Creates function to build the model that we'll run in VRPSolver.
function build_model(data::DataOPHS, app)
    total_cuts = 0 #Initializes the total number of cuts used.
    E = edges(data) # Set of edges of the input graph G′. It is the same of data.G'.E
    H = data.Hotels # Set of hotels. The tour starts at H[1] and ends at hotel H[2].
    C = data.Customers # Set of customers.
    n = length(H)+length(C) # Total number of vertices.
    Tmax = data.Tmax # Total tour length.
    D = data.D # Number of trips.
    Td = data.Td # Vector containing the length of each trip. If data.symmetric == true then the maximum length of each trip is equal to Td[1]
    
    T=[i for i=1:D] # Set of trip indexes.
    
    # Function ed will define which direction will be mapped 
    ed(i, j) = i < j ? (i, j) : (j, i)

    
    #= Lets start building the model. 
        Case 1 -> data is symmetric (the maximum length is the same in the trips)
        Case 2 -> data is not symmetric.
        Case 3 -> A single graph with stops
    =#

    # Formulation
    bwtsp = VrpModel()

    # Variables
    
    @variable(bwtsp.formulation, x[e in E], Int) # Number of times edge e is traversed.
    @variable(bwtsp.formulation, 0 <= g[i in C] <= 1, Int) # g_i = 1 if the customer i is visited.
    
    # If symmetric, we don't need to keep the moment a trip was made. The order does not influence the result.
    if data.symmetric
        @variable(bwtsp.formulation, y[i in H], Int) # y[i] is the number of times that the hotel i is visited
    else
        @variable(bwtsp.formulation, y[i in H, k in T], Int) # y[i,k] = 1 if trip k starts at the hotel i
        @variable(bwtsp.formulation, t[i in H, k in T], Int) # t[i,k] = 1 if trip k ends at the hotel i
    end
    
    # Objective Function
    # The object function maximizes the total score obtained with the visited customers
    @objective(bwtsp.formulation, Min, -sum(s(data,i)*g[i] for i in C))

    # Constraints
    
    # The total tour length does not exceed Tmax
    if data.symmetric
        # If symmetric, we sum the cost of a "dummy arc" from the start to the end
        # This arc was created for modeling purposes. It does not exist in fact. 
        @constraint(bwtsp.formulation, sum(c(data,e)*x[e] for e in E) <= Tmax + c(data,(1,2))) 
    else 
        if !app["complete"]
            @constraint(bwtsp.formulation, sum(c(data,e)*x[e] for e in E) <= Tmax)
        end 
    end

    # Is each customer visited only once? 
    # Is the degree of a visited customer 2?
    @constraint(bwtsp.formulation, deg[i in C], sum(x[e] for e in E if e[1] == i || e[2]==i) == 2*g[i]) 
    
    @constraint(bwtsp.formulation, men[i in C], g[i] <= 1) 
    @constraint(bwtsp.formulation, mai[i in C], g[i] >= 0)
    
    
    if data.symmetric
        @constraint(bwtsp.formulation, y[H[1]] >= 1) # H[1] is visited at least one time
        @constraint(bwtsp.formulation, y[H[2]] >= 1) # H[2] is visited at least one time
        @constraint(bwtsp.formulation, x[(1,2)] >= 1) # Edge (1,2) is always traversed at least once because we consider the dummy trip H[2] -> H[1].
        @constraint(bwtsp.formulation,[i in H],sum(x[e] for e in δ(data, i) if e in E) == 2*y[i]) # degree of the hotel
    else
        @constraint(bwtsp.formulation, y[H[1],1] == 1) #The first trip starts at H[1]
        @constraint(bwtsp.formulation, t[H[2],D] == 1) #The last trip ends at H[2]
        #@constraint(bwtsp.formulation, hotel_deg[i in H], sum(x[e] for e in δ(data, i) if e in E) == sum(y[i,k] for k in T) + sum(t[i,k] for k in T)) #Degree of each hotel
        @constraint(bwtsp.formulation, degt[k = 2:length(T) , i in H], y[i,k] == t[i,k-1]) # The k-th trip starts at the same hotel where the k-1 th trip ends
    end
    
    #println(bwtsp.formulation)
    # Now we'll build the graphs.
    
    # A small that could be used to represent the dummy trip
    # This is not used currently
    function build_trip_1_2_graph()
        V′ = [1,2] # just the two vertices are necessary. This graph is only for the dummy trip
        v_source, v_sink = 1, 2
        L, U = 1,1 # The graph is used to generate a single trip    
        G = VrpGraph(bwtsp, V′, v_source, v_sink, (L, U))
        arc_id = add_arc!(G, v_source, v_sink)
        add_arc_var_mapping!(G, arc_id, [x[(1,2)],y[1],y[2]])
        return G
    end
    

    # The non_symmetric case.
    
    #=
    In the non_symmetric case, we build a graph for each trip.
    As the trips have different lengths, we give it as input. 
    =#

    function build_graph_non_symmetric(trip_index, trip_length)
        #= 
        # The Vertices of the graph are:
            # {0} ⋃ Customers ⋃ Hotels ⋃ copies of Hotels ⋃ {|C|+|H|+|H|+1}
            # where:
                # 0 is the v_source
                # |C|+|H|+|H|+1 is the v_sink  
        =#
        
        V′ = [i for i in 0:length(C)+length(H)+length(H)+1]   
        v_source, v_sink = 0, length(C)+length(H)+length(H)+1
        L, U = 1,1 # Each graph is used to generate a single trip
        G = VrpGraph(bwtsp, V′, v_source, v_sink, (L, U))

        # Resouces
        # Resource used to ensure the trip length constraint
        dist_res_id = add_resource!(G, main=true, disposable = true) 
        
        # Setting the resource bounds
        for i in V′
            # All vertices have the same bounds, between 0 and the trip_length 
            set_resource_bounds!(G, i, dist_res_id, 0.0, trip_length)
        end

        # Mapping arcs and setting arcs` consumption

        # Arcs from v_source to each hotel and from each hotel`s copy to the v_sink
        for i in H
            # from the v_source to the hotel i
            arc_id = add_arc!(G, v_source, i)
            add_arc_var_mapping!(G, arc_id, y[i,trip_index])
            set_arc_consumption!(G, arc_id, dist_res_id, 0.0)
            
            # from the copy of hotel i to v_sink
            arc_id = add_arc!(G, i+n, v_sink)
            add_arc_var_mapping!(G, arc_id, t[i,trip_index])
            set_arc_consumption!(G, arc_id, dist_res_id, 0.0)
        end
        
        for (i,j) in E
            #arcs between customers
            if i in C && j in C
                arc_id = add_arc!(G, i, j)
                add_arc_var_mapping!(G, arc_id, x[(i, j)])
                set_arc_consumption!(G, arc_id, dist_res_id, d(data,(i,j)))
    
                arc_id = add_arc!(G, j, i)
                add_arc_var_mapping!(G, arc_id, x[(i, j)])
                set_arc_consumption!(G, arc_id, dist_res_id, d(data,(i,j)))
            end

            # Arcs between a hotel and a customer and between a customer and a hotel copy
            if i in H && j in C
                #from the hotel i to the customer j
                arc_id = add_arc!(G, i, j)
                add_arc_var_mapping!(G, arc_id, x[(i, j)])
                set_arc_consumption!(G, arc_id, dist_res_id, d(data,(i,j)))
    
                # from the the customer j and the copy of hotel i, given by the index i+n
                arc_id = add_arc!(G, j, i+n)
                add_arc_var_mapping!(G, arc_id, x[(i, j)])
                set_arc_consumption!(G, arc_id, dist_res_id, d(data,(i,j)))
            end
            
            #arcs between hotels and copies
            if i in H && j in H
                # from a hotel i to the copy of a hotel j (given by j+n)
                arc_id = add_arc!(G, i, j+n)
                add_arc_var_mapping!(G, arc_id, x[ed(i, j)])
                set_arc_consumption!(G, arc_id, dist_res_id, d(data,ed(i,j)))

                # from a hotel j to the copy of a hotel i (given by i+n)
                # Parece que estamos mapeando de novo por cima...

                arc_id = add_arc!(G, j, i+n)
                add_arc_var_mapping!(G, arc_id, x[ed(i, j)])
                set_arc_consumption!(G, arc_id, dist_res_id, d(data,ed(i,j)))
            end
        end
        return G
    end
    
    #= The symmetric case.
        Now, we dont need the input trip_length as it is the same for all the trips
        All the trips will occur in the smae graph.
    =# 

    function build_graph_symmetric()
        #= 
        # The Vertices of the graph are:
            # {0} ⋃ Customers ⋃ Hotels ⋃ copies of Hotels ⋃ {|C|+|H|+|H|+1}
            # where:
                # 0 is the v_source
                # |C|+|H|+|H|+1 is the v_sink  
        =#
        V′ = [i for i in 0:length(C)+length(H)+length(H)+1] # {0} ⋃ Customers ⋃ Hotels ⋃ copies of Hotels ⋃ {|C|+|H|+|H|+1}, where the vertices 0 and |C|+|H|+|H|+1 are dummy vertices used to indicate repsctively the v_source and the v_sink  
        v_source, v_sink = 0, length(C)+length(H)+length(H)+1
        L, U = D+1,D+1 # Adding 1 to the number of trips because of the dummy trip H[2]->H[1]  
        G = VrpGraph(bwtsp, V′, v_source, v_sink, (L, U))
        
        # Resouces
        # Resource used to ensure the trip length constraint
        dist_res_id = add_resource!(G, main=true)
        
        for i in V′ 
            # All the vertices have the same bounds. Now, the upper limit is Td[1], which is the same of Td[2]...
            set_resource_bounds!(G, i, dist_res_id, 0.0, Td[1])
        end

        # Mapping arcs and setting arcs` consumption
        # Arcs from v_source to each hotel and from each hotel`s copy to the v_sink
        for i in H
            #from the v_source to the hotel i
            arc_id = add_arc!(G, v_source, i)
            add_arc_var_mapping!(G, arc_id, y[i])
            set_arc_consumption!(G, arc_id, dist_res_id, 0.0)
            
            #from the hotel i copy to the v_sink
            arc_id = add_arc!(G, i+n, v_sink)
            set_arc_consumption!(G, arc_id, dist_res_id, 0.0)
        end

        for (i,j) in E
            if i != 0 && j != 0
                #arcs between customers
                if i in C && j in C
                    arc_id = add_arc!(G, i, j)
                    add_arc_var_mapping!(G, arc_id, [x[(i, j)]])
                    set_arc_consumption!(G, arc_id, dist_res_id, d(data,(i,j)))
        
                    arc_id = add_arc!(G, j, i)
                    add_arc_var_mapping!(G, arc_id, [x[(i, j)]])
                    set_arc_consumption!(G, arc_id, dist_res_id, d(data,(i,j)))
                end
                
                #arcs between a hotel and a customer and between a customer and a hotel copy
                if i in H && j in C
                    arc_id = add_arc!(G, i, j)
                    add_arc_var_mapping!(G, arc_id, [x[(i, j)]])
                    set_arc_consumption!(G, arc_id, dist_res_id, d(data,(i,j)))
        
                    arc_id = add_arc!(G, j, i+n)
                    add_arc_var_mapping!(G, arc_id, x[(i, j)])
                    set_arc_consumption!(G, arc_id, dist_res_id, d(data,(i,j)))
                end
                
                if i in H && j in H
                    #from a hotel i to the copy of a hotel j
                    arc_id = add_arc!(G, i, j+n)
                    add_arc_var_mapping!(G, arc_id, x[ed(i, j)])
                    set_arc_consumption!(G, arc_id, dist_res_id, d(data,ed(i,j)))
                    
                    #from a hotel j to the copy of a hotel i
                    arc_id = add_arc!(G, j, i+n)
                    add_arc_var_mapping!(G, arc_id, x[ed(i, j)])
                    set_arc_consumption!(G, arc_id, dist_res_id, d(data,ed(i,j)))
                end                
            end
        end
        return G
    end

     # The alternative graph for the symmetric case.
    function build_alternative_graph_symmetric()
        #= 
        # The Vertices of the graph are:
            # {0} ⋃ Customers ⋃ Hotels ⋃ copies of Hotels ⋃ {|C|+|H|+|H|+1}
            # where:
                # 0 is the v_source
                # |C|+|H|+|H|+1 is the v_sink  
        =#

        V′ = vcat([0],C) # {0} ⋃ Customers
        v_source, v_sink = 0, 0
        
        L, U = D,D # Adding 1 to the number of trips because of the dummy trip H[2]->H[1]  
        G = VrpGraph(bwtsp, V′, v_source, v_sink, (L, U))
        
        # Resouces
        # Resource used to ensure the trip length constraint
        dist_res_id = add_resource!(G, main=true) #Resource used to ensure the trip length constraint
        for i in V′
            l_i, u_i = 0.0, Td[1] 
            set_resource_bounds!(G, i, dist_res_id, l_i, u_i)
        end



        # Mapping arcs and setting arcs` consumption
        
        # Arcs from v_source to each hotel and from each hotel`s copy to the v_sink
        for j in C
            #from the v_source to the hotel i
            for i in H
                arc_id = add_arc!(G, v_source, j)
                add_arc_var_mapping!(G, arc_id, [x[(i, j)]])
                set_arc_consumption!(G, arc_id, dist_res_id, d(data,(i,j)))

                arc_id = add_arc!(G, j, v_sink)
                add_arc_var_mapping!(G, arc_id, [x[(i, j)]])
                set_arc_consumption!(G, arc_id, dist_res_id, d(data,(i,j)))
            end
        end

        for i in C
            for j in C
                if i < j
                    arc_id = add_arc!(G, i, j)
                    add_arc_var_mapping!(G, arc_id, [x[(i, j)]])
                    set_arc_consumption!(G, arc_id, dist_res_id, d(data,(i,j)))

                    arc_id = add_arc!(G, j, i)
                    add_arc_var_mapping!(G, arc_id, [x[(i, j)]])
                    set_arc_consumption!(G, arc_id, dist_res_id, d(data,(i,j)))
                end
            end
        end

        return G
    end    

    # The symmetric case.
    function build_graph_CompleteRoute_symmetric()
        #= 
        # The Vertices of the graph are:
            # Customers ⋃ Hotels
        =#       
        V′ = [i for i in 1:length(C)+length(H)]   
        v_source, v_sink = 1, 2
        L, U = 1,1 # Adding 1 to the number of trips because of the dummy trip H[2]->H[1]  
        G = VrpGraph(bwtsp, V′, v_source, v_sink, (L, U))
        
        # Resouces
        # Resource used to ensure the trip length constraint
        time_res_id = add_resource!(G, main=true)
        dist_res_id = add_resource!(G, main=true) #Resource used to ensure the trip length constraint
        
        for i in V′
            i in H ? set_resource_bounds!(G, i, dist_res_id, 0.0, 0.0) : set_resource_bounds!(G, i, dist_res_id, 0.0, Td[1])
            set_resource_bounds!(G, i, time_res_id, 0.0, Tmax)
        end
        
        # Mapping arcs and setting arcs` consumption
        
        # Arcs from hotels to  to each hotel and from each hotel to v_sink
        
        for e in E
            arc_id = add_arc!(G, e[1], e[2])
            set_arc_consumption!(G, arc_id, time_res_id, d(data,e))
            if e[2] in H
                add_arc_var_mapping!(G, arc_id, [x[e], y[e[2]]])                
                set_arc_consumption!(G, arc_id, dist_res_id, d(data,e)-Td[1])
            else
                add_arc_var_mapping!(G, arc_id, [x[e]])
                set_arc_consumption!(G, arc_id, dist_res_id, d(data,e))
            end

            arc_id = add_arc!(G, e[2], e[1])
            set_arc_consumption!(G, arc_id, time_res_id, d(data,e))
            if e[1] in H
                add_arc_var_mapping!(G, arc_id, [x[e], y[e[1]]])                
                set_arc_consumption!(G, arc_id, dist_res_id, d(data,e)-Td[1])
            else
                add_arc_var_mapping!(G, arc_id, [x[e]])
                set_arc_consumption!(G, arc_id, dist_res_id, d(data,e))
            end

        end    
        return G
    end

    #=
    # The case considering stops with negative consumption.
    function build_graph_stops()
        #= 
        # For each trip, we have 2 hotels' layers and one customer's layer.
            # {0} ⋃ Customers * D ⋃ Hotels * D + 1 ⋃ {|C|*D + |H|*D+1 + 1} 
        =#
        
        C_layers =[]
        H_layers =[]
        # i is the layer
        for i = 1:D
            first_id = length(C) * (i-1)  + 1
            last_id = length(C) * i
            push!(C_layers, Layer("Customer", i, first_id, last_id))
            
            first_id = D * length(C) + length(H) * (i-1)  + 1
            last_id = D * length(C) + length(H) * i 
            push!(H_layers, Layer("Hotel", i, first_id, last_id))
            
            # an extra hotel l
            ayer 
            if i == D
                first_id = D * length(C) + length(H) * (D)  + 1
                last_id = D * length(C) + length(H) * (D+1)
                push!(H_layers, Layer("Hotel", D+1, first_id, last_id))
            end
        end


        V′ = [i for i in 0:H_layers[end].last_id]  
        v_source, v_sink = 0, H_layers[end].last_id
        
        L, U = D,D  
        G = VrpGraph(bwtsp, V′, v_source, v_sink, (L, U))
        
        # Defining the resources
        
        # Resource used to ensure the trip length constraint
        dist_res_id = add_resource!(G, main = true, disposable = true) 
        
        # Now we have use binary resources to ensure the lenght trips are ok
        bin_res_id = []
        for i = 1:D
            push(bin_res_id, add_resource!(G, main = false, binary = true))
        end
        
        
        
        for i in V′
            l_i, u_i = 0.0, Td[1] 
            set_resource_bounds!(G, i, dist_res_id, l_i, u_i)
        end

        # Arcs from the v_source to each hotel and from each copy of hotel to the v_sink
        for i in H
            #from the v_source to the hotel i
            arc_id = add_arc!(G, v_source, i)
            add_arc_var_mapping!(G, arc_id, y[i])
            set_arc_consumption!(G, arc_id, dist_res_id, 0.0)
            
            #from the hotel i copy to the v_sink
            arc_id = add_arc!(G, i+n, v_sink)
            set_arc_consumption!(G, arc_id, dist_res_id, 0.0)
        end

        for (i,j) in E
            if i != 0 && j != 0
                #arcs between customers
                if i in C && j in C
                    arc_id = add_arc!(G, i, j)
                    add_arc_var_mapping!(G, arc_id, [x[(i, j)]])
                    set_arc_consumption!(G, arc_id, dist_res_id, d(data,(i,j)))
        
                    arc_id = add_arc!(G, j, i)
                    add_arc_var_mapping!(G, arc_id, [x[(i, j)]])
                    set_arc_consumption!(G, arc_id, dist_res_id, d(data,(i,j)))
                end
                
                #arcs between a hotel and a customer and between a customer and a hotel copy
                if i in H && j in C
                    arc_id = add_arc!(G, i, j)
                    add_arc_var_mapping!(G, arc_id, [x[(i, j)]])
                    set_arc_consumption!(G, arc_id, dist_res_id, d(data,(i,j)))
        
                    arc_id = add_arc!(G, j, i+n)
                    add_arc_var_mapping!(G, arc_id, x[(i, j)])
                    set_arc_consumption!(G, arc_id, dist_res_id, d(data,(i,j)))
                end
                
                if i in H && j in H
                    #from a hotel i to the copy of a hotel j
                    arc_id = add_arc!(G, i, j+n)
                    add_arc_var_mapping!(G, arc_id, x[ed(i, j)])
                    set_arc_consumption!(G, arc_id, dist_res_id, d(data,ed(i,j)))
                    
                    #from a hotel j to the copy of a hotel i
                    arc_id = add_arc!(G, j, i+n)
                    add_arc_var_mapping!(G, arc_id, x[ed(i, j)])
                    set_arc_consumption!(G, arc_id, dist_res_id, d(data,ed(i,j)))
                end                
            end
        end
        return G
    end
    =#

    graphs =[]
    if data.symmetric
        if app["complete"]
            G = build_graph_CompleteRoute_symmetric()
        else
            #G_= build_trip_1_2_graph()
            #add_graph!(bwtsp, G_)
            G = build_graph_symmetric()
            #G = build_alternative_graph_symmetric()
        end
        add_graph!(bwtsp, G)
        push!(graphs, G)
    
    else
        for k=1:length(T)
            if k != length(T)+1
                G = build_graph_non_symmetric(k,Td[k]+0.01) 
                push!(graphs, G)
            else
                G = build_trip_1_2_graph()
            end
            add_graph!(bwtsp, G)            
        end
    end

    set_vertex_packing_sets!(bwtsp, [[(G,i) for G in graphs] for i in C])

    for G in graphs
        define_elementarity_sets_distance_matrix!(bwtsp, G, [[d(data,ed(i, j)) for i in C] for j in C])
    end
    
    set_branching_priority!(bwtsp, "x", 1)
    set_branching_priority!(bwtsp, "y", 1)
    if !data.symmetric
        set_branching_priority!(bwtsp, "t", 1)
    end
    set_branching_priority!(bwtsp, "g", 1)

    function maxflow_mincut_callback()
        M = 100000
        graph = SparseMaxFlowMinCut.ArcFlow[]
        especial_flow = 0.0
        for (i, j) in E
            e = (i, j)
            value::Float64 = get_value(bwtsp.optimizer, x[e])
            
            if value > 0.0001
                flow_::Int = trunc(floor(value, digits=5) * M)
                push!(graph, SparseMaxFlowMinCut.ArcFlow(i, j, flow_)) # arc i -> j
                push!(graph, SparseMaxFlowMinCut.ArcFlow(j, i, flow_)) # arc j -> i
            end
        end
        
        added_cuts = []
        s = H[1]
        for c in C
            value::Float64 = get_value(bwtsp.optimizer, g[c])
            if value > 0.001
                maxFlow, flows, cut = SparseMaxFlowMinCut.find_maxflow_mincut(SparseMaxFlowMinCut.Graph(n, graph), s, c)
                if (maxFlow / M) < (2*value - 0.001) && !in(cut, added_cuts)
                    set1, set2 = [], []
                    [cut[i] == 1 ? push!(set1, i) : push!(set2, i) for i in 1:n]
                    #@show set1, set2
                    
                    # A is an array of variables that will be in a cut
                    A = [x[ed(i, j)] for i in set1 for j in set2 if ed(i, j) in E]
                    push!(A, g[c])
                    # B is an array of coefficients
                    B = [1.0 for i in set1 for j in set2 if ed(i, j) in E]
                    push!(B,-2.0)
                    add_dynamic_constr!(bwtsp.optimizer, A, B, >=, 0.0, "mincut")
                    push!(added_cuts, cut)
                end
            end
        end
        if length(added_cuts) > 0
            data.cuts += length(added_cuts)
            @show data.cuts
            println(">>>>> Add min cuts : ", length(added_cuts), " cut(s) added") 
        end
    end

    if data.symmetric && !app["complete"]
        add_cut_callback!(bwtsp, maxflow_mincut_callback, "mincut")
    end
 
    if data.symmetric
        return (bwtsp, x,y,g)
    else
        return (bwtsp, x,y,t,g)
    end
end