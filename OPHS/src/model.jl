function build_model(data::DataOPHS)
    total_cuts = 0 #Initializes the total number of cuts 
    E = edges(data) #Set of edges of the input graph G′
    H = data.Hotels #Set of hotels. The tour starts at H[1] and ends at hotel H[2]
    C = data.Customers #Set of customers
    n = length(H)+length(C) #Total number of vertices
    Tmax = data.Tmax #Total tour length
    D = data.D #Number of trips
    Td = data.Td # Vector containing the length of each trip. If data.symmetric == true then the maximum length of each trip is equal to Td[1]
    T=[i for i=1:D] # Set of trip indexes.
    
    ed(i, j) = i < j ? (i, j) : (j, i)

    # Formulation
    bwtsp = VrpModel()

    # Variables
    @variable(bwtsp.formulation, x[e in E],Int) # Number of times that a given edge e is traversed 
    @variable(bwtsp.formulation, 0 <= g[i in C] <= 1, Int) # g_i = 1 if the customer i is visited
    if data.symmetric
        @variable(bwtsp.formulation, y[i in H], Int) #Number of times that the hotel i is visited
    else
        @variable(bwtsp.formulation, y[i in H, k in T], Int) #y[i,k] = 1 if the trip k starts at the hotel i
        @variable(bwtsp.formulation, t[i in H, k in T], Int) #t[i,k] = 1 if the trip k ends at the hotel i
    end
    
    # Objective Function
    @objective(bwtsp.formulation, Min, -sum(s(data,i)*g[i] for i in C))

    # Constraints
    if data.symmetric
        @constraint(bwtsp.formulation, sum(c(data,e)*x[e] for e in E) <= Tmax + c(data,(1,2))) # The total tour length does not exceed Tmax
    else 
        @constraint(bwtsp.formulation, sum(c(data,e)*x[e] for e in E) <= Tmax) # The total tour length does not exceed Tmax
    end

    @constraint(bwtsp.formulation, deg[i in C], sum(x[e] for e in E if e[1] == i || e[2]==i) == 2*g[i]) # degree of each customer
    @constraint(bwtsp.formulation, men[i in C],g[i] <= 1)
    @constraint(bwtsp.formulation, mai[i in C],g[i] >= 0)
    if !data.symmetric
        @constraint(bwtsp.formulation, y[H[1],1] == 1) #The first trip starts at H[1]
        @constraint(bwtsp.formulation, t[H[2],D] == 1) #The last trip ends at H[2]
        #@constraint(bwtsp.formulation, hotel_deg[i in H], sum(x[e] for e in δ(data, i) if e in E) == sum(y[i,k] for k in T) + sum(t[i,k] for k in T)) #Degree of each hotel
        @constraint(bwtsp.formulation, degt[k=2:length(T),i in H], y[i,k] == t[i,k-1]) # The k-th trip starts at the same hotel where the k-1 th trip ends 
    else
        @constraint(bwtsp.formulation, y[H[1]] >= 1) # H[1] is visited at least one time
        @constraint(bwtsp.formulation, y[H[2]] >= 1) # H[2] is visited at least one time
        @constraint(bwtsp.formulation, x[(1,2)] >= 1) # the edge (1,2) is always traversed at least one time because we consider the dummy trip H[2] -> H[1].
        @constraint(bwtsp.formulation,[i in H],sum(x[e] for e in δ(data, i) if e in E) == 2*y[i]) # degree of the hotel
    end
    
    #println(bwtsp.formulation)

    # Build the model directed graph G=(V,A)
    function build_graph_non_symmetric(trip_index, trip_length)
        V′ = [i for i in 0:length(C)+length(H)+length(H)+1] # {0} ⋃ Customers ⋃ Hotels ⋃ copies of Hotels ⋃ {|C|+|H|+|H|+1}, where the vertices 0 and |C|+|H|+|H|+1 are dummy vertices used to indicate repsctively the v_source and the v_sink  
        v_source, v_sink = 0, length(C)+length(H)+length(H)+1
        
        L, U = 1,1 # Each graph is used to generate a single trip

        G = VrpGraph(bwtsp, V′, v_source, v_sink, (L, U))
        dist_res_id = add_resource!(G, main=true) #Resource used to ensure the trip length constraint
        
        for i in V′
            l_i, u_i = 0.0, trip_length
            set_resource_bounds!(G, i, dist_res_id, l_i, u_i)
        end

        # Arcs from the v_source to each hotel and from each copy of hotel to the v_sink
        for i in H
            #from the v_source to the hotel i
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

            #arcs between a hotel and a customer and between a customer and a hotel copy
            if i in H && j in C
                #from the hotel i to the customer j
                arc_id = add_arc!(G, i, j)
                add_arc_var_mapping!(G, arc_id, x[(i, j)])
                set_arc_consumption!(G, arc_id, dist_res_id, d(data,(i,j)))
    
                #from the the customer j and the copy of hotel i
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
        return G
    end

    function build_graph_symmetric()
        V′ = [i for i in 0:length(C)+length(H)+length(H)+1] # {0} ⋃ Customers ⋃ Hotels ⋃ copies of Hotels ⋃ {|C|+|H|+|H|+1}, where the vertices 0 and |C|+|H|+|H|+1 are dummy vertices used to indicate repsctively the v_source and the v_sink  
        v_source, v_sink = 0, length(C)+length(H)+length(H)+1
        
        L, U = D+1,D+1 # We are increasing 1 to the number of trips because we consider the dummy trip H[2]->H[1]  

        G = VrpGraph(bwtsp, V′, v_source, v_sink, (L, U))
        dist_res_id = add_resource!(G, main=true) #Resource used to ensure the trip length constraint
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

    graphs =[]
    if data.symmetric
        G = build_graph_symmetric()
        add_graph!(bwtsp, G)
        push!(graphs, G)
    else
        for k in T
            G = build_graph_non_symmetric(k,Td[k]+0.01) 
            add_graph!(bwtsp, G)
            push!(graphs, G)
        end
    end

    set_vertex_packing_sets!(bwtsp, [[(G,i) for G in graphs] for i in C])

    for G in graphs
        define_elementarity_sets_distance_matrix!(bwtsp, G, [[d(data,ed(i, j)) for i in C] for j in C])
    end
    
    set_branching_priority!(bwtsp, "x", 1)
    set_branching_priority!(bwtsp, "y", 2)
    if !data.symmetric
        set_branching_priority!(bwtsp, "t", 2)
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
                    A = [x[ed(i, j)] for i in set1 for j in set2 if ed(i, j) in E]
                    push!(A, g[c])
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

    if data.symmetric
        add_cut_callback!(bwtsp, maxflow_mincut_callback, "mincut")
    end
 
    if data.symmetric
        return (bwtsp, x,y,g)
    else
        return (bwtsp, x,y,t,g)
    end
end