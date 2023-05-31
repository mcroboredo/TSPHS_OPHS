import Unicode

mutable struct Vertex
   id_vertex::Int
   pos_x::Float64
   pos_y::Float64
   s_time::Float64 # service time
end

# Undirected graph
mutable struct InputGraph
   V′::Array{Vertex} # set of vertices
   E::Array{Tuple{Int64,Int64}} # set of edges
   cost::Dict{Tuple{Int64,Int64},Float64} # cost for each edge
end

mutable struct DataTSPHS
   G′::InputGraph # input graph
   H′::Array{Int64} # hotels set
   C′::Array{Int64} # customers set
   Lim::Int # time limit
end

function distance(data::DataTSPHS, arc::Tuple{Int64, Int64})
   e = (arc[1] < arc[2]) ? arc : (arc[2],arc[1])
   if haskey(data.G′.cost, e) # use already calculated value
      return data.G′.cost[e]
   elseif e[1] == e[2]
      return 0.0
   else
      u, v = arc
      vertices = data.G′.V′
      # euclidian distance
      x_sq = (vertices[v].pos_x - vertices[u].pos_x)^2
      y_sq = (vertices[v].pos_y - vertices[u].pos_y)^2
      return floor(sqrt(x_sq + y_sq), digits=1)
   end
end



function readTSPHSData(app::Dict{String,Any})

   str = Unicode.normalize(read(app["instance"], String); stripcc=true)
   breaks_in = [' '; ':'; '\n']
   aux = split(str, breaks_in; limit=0, keepempty=false)

   G′ = InputGraph([],[],Dict())
   data = DataTSPHS(G′, [], [], 0)

   nodes = []
   h, c, l = 0, 0, 0

   for i in 1:length(aux)
      h = parse(Int, aux[i])
      c = parse(Int, aux[i+1])
      l = parse(Int, aux[i+2])

      j = i + 3
      last = 1

      while last <= h # hotels vertices
         v = Vertex(0, 0, 0, 0)
         v.id_vertex = last
         v.pos_x = parse(Float64, aux[j+1])
         v.pos_y = parse(Float64, aux[j+2])
         v.s_time = 0.0
         push!(nodes, v) # add v in the vertex array
         last+=1
         j+=3
      end

      while last <= (h + c) # customers vertices
         v = Vertex(0, 0, 0, 0)
         v.id_vertex = last
         v.pos_x = parse(Float64, aux[j+1])
         v.pos_y = parse(Float64, aux[j+2])
         v.s_time = parse(Float64, aux[j+3])
         push!(nodes, v) # add v in the vertex array
         last+=1
         j+=4
      end

      if j > length(aux)
         break
      end
   end

   data.H′ = [i.id_vertex for i in nodes[1:h]]
   data.C′ = [i.id_vertex for i in nodes[h + 1:end]]
   data.Lim = l
   G′.V′ = nodes # add vertices to graph G′
   # println(G′.V′)

   for i in vertices(data)
      for j in vertices(data) # add edges between customers and hotels
         if i < j
            e = (i,j)
            push!(G′.E, e) # add edge e
            data.G′.cost[e] = distance(data, e) # cost edge e
         end
      end
   end
   # println(data.G′.V′)
   # println(data.G′.cost)

   return data
end

edges(data::DataTSPHS) = data.G′.E # return set of edges
c(data,e) = data.G′.cost[e] # cost of the edge e
service(data,i) = data.G′.V′[i].s_time # service time of vertex i
dimension(data::DataTSPHS) = length(data.G′.V′) # return number of vertices
vertices(data::DataTSPHS) = [i.id_vertex for i in data.G′.V′[1:end]] # return set of vertices (hotels+customers)
nb_customers(data::DataTSPHS) = length(data.C′) # return number of customers
nb_hotels(data::DataTSPHS) = length(data.H′) # return number of hotels

# return incident edges of i
function δ(data::DataTSPHS, i::Integer)
   incident_edges = Vector{Tuple}()
   for j in 1:i-1 push!(incident_edges, (j, i)) end
   for j in i+1:(length(data.G′.V′)) push!(incident_edges, (i, j)) end
   return incident_edges
end
