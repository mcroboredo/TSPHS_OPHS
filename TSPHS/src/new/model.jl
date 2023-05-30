# using SimpleWeightedGraphs
function build_model(data::DataTSPHS, app::Dict{String,Any})

   E′ = edges(data) # set of edges of the input graph G′
   n = dimension(data) # number of vertices
   H = data.H′ # Set of hotels
   C = data.C′ # Set of customers
   q = app["qvalue"] # q value

   # set of customers of the input graph G′ + source and sink (0)
   V = [i for i in 0:n]

   #v(i) = (i <= n) ? i : i - n # convert hotel vertex duplicates to original hotel vertex
   ed(i,j) = i < j ? (i,j) : (j,i) # get the proper edge
   d((i,j),data) = i != j ? c(data,(i,j)) : 0.0 # get the distance between i and j

   # Formulation
   tsphs = VrpModel()
   @variable(tsphs.formulation, x[e in E′], Int)
   @variable(tsphs.formulation, 0 <= b[h in H] <= q, Int)
   @objective(tsphs.formulation, Min, sum(c(data,e) * x[e] for e in E′))
   @constraint(tsphs.formulation, degC[i in C], sum(x[e] for e in δ(data, i)) == 2.0)
   @constraint(tsphs.formulation, degH[h in H], sum(x[e] for e in δ(data, h)) == 2.0 * b[h])
   @constraint(tsphs.formulation, b[1] >= 1)
   # @constraint(tsphs.formulation, sum(b[h] for h in H) == q)
   # println(tsphs.formulation)

   # Build the model directed graph G=(V,A)
   function build_graph()

      v_source = v_sink = 0
      L = U = q

      # node ids of G from 0 to |V|
      G = VrpGraph(tsphs, V, v_source, v_sink, (L, U))
      cap_res_id = add_resource!(G, main = true) # R = R_M = {cap_res_id}
      for i in V
         l_i, u_i = 0.0, Float64(data.Lim) # accumulated resource consumption interval [l_i, u_i] for the vertex i
         set_resource_bounds!(G, i, cap_res_id, l_i, u_i)
      end

      # Build set of arcs A from E′ (two arcs for each edge (i,j))
      for i in H # setting the arcs between source, sink, and hotels
         arc_id = add_arc!(G, v_source, i) # source -> i(Hotel)
         set_arc_consumption!(G, arc_id, cap_res_id, 0.0)
         arc_id = add_arc!(G, i, v_sink) # i(Hotel) -> sink
         set_arc_consumption!(G, arc_id, cap_res_id, 0.0)
      end

      for (i,j) in E′
         service_cost = (service(data,i) + service(data,j))/2
         c_e = c(data,(i,j)) + service_cost

         # add arcs i - > j
         arc_id = add_arc!(G, i, j)
         add_arc_var_mapping!(G, arc_id, x[(i,j)])
         set_arc_consumption!(G, arc_id, cap_res_id, c_e)
         # add arcs j - > i
         arc_id = add_arc!(G, j, i)
         add_arc_var_mapping!(G, arc_id, x[(i,j)])
         set_arc_consumption!(G, arc_id, cap_res_id, c_e)
      end

      return G
   end

   G = build_graph()
   add_graph!(tsphs, G)
   # println(G)

   set_vertex_packing_sets!(tsphs, [[(G,i)] for i in C])

   define_elementarity_sets_distance_matrix!(tsphs, G, [[ distance(data, (i,j)) for i in C] for j in C])

   set_branching_priority!(tsphs, "b", 1)
   set_branching_priority!(tsphs, "x", 1)

   function maxflow_mincut_callback()
      M = 100000
      g = SparseMaxFlowMinCut.ArcFlow[]
      for (i,j) in E′
         e = (i,j)
         value::Float64 = get_value(tsphs.optimizer, x[e])
         if  value > 0.0001
            flow_::Int = trunc(floor(value, digits=5) * M)
            push!(g, SparseMaxFlowMinCut.ArcFlow(i, j, flow_)) # arc i -> j
            push!(g, SparseMaxFlowMinCut.ArcFlow(j, i, flow_)) # arc j -> i
         end
      end

      added_cuts = []
      s = H[1]
      for t in 1:length(C)
         maxFlow, flows, cut = SparseMaxFlowMinCut.find_maxflow_mincut(SparseMaxFlowMinCut.Graph(n, g), s, C[t])
         if (maxFlow/M) < (2 - 0.001) && !in(cut, added_cuts)
            set1, set2 = [], []
            [cut[i] == 1 ? push!(set1, i) : push!(set2, i) for i in 1:n]
            add_dynamic_constr!(tsphs.optimizer, [x[ed(i,j)] for i in set1 for j in set2], [1.0 for i in set1 for j in set2], >=, 2.0, "mincut")
            push!(added_cuts, cut)
         end
      end
      length(added_cuts) > 0 && printstyled(">>>>> Add min cuts : $(length(added_cuts)) cut(s) added\n", color=:yellow)
   end
   add_cut_callback!(tsphs, maxflow_mincut_callback, "mincut")

   return (tsphs, x)
end
