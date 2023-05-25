import Unicode

# Struct of a Vertex (customer and hotel)  
mutable struct Vertex
    id_vertex::Int # id of the vertex
    pos_x::Float64 # position at x coord
    pos_y::Float64 # position at y coord
    score::Float64 # Score of each customer
end

# Undirected graph
mutable struct InputGraph
    V′::Array{Vertex} # set of vertices (Customers U Hotels). Each vertex is a struct as above
    E::Array{Tuple{Int64,Int64}} # set of edges
    cost::Dict{Tuple{Int64,Int64},Float64} # cost for each edge
end

mutable struct DataOPHS
    G′::InputGraph
    Hotels::Array{Int64} # Hotel nodes. The first two hotels are respectively the start depot and the end depot 
    Customers::Array{Int64} # Customer nodes
    Tmax::Float64 # Total tour length
    D::Int64 # Number of trips
    Td::Array{Float64} # array containing the length of each trip
    symmetric::Bool # = true if Td[1] = Td[2] = ... = Td[D]
    cuts::Int64 #Total number of cuts inserted by the user callback
end

vertices(data::DataOPHS) = [i.id_vertex for i in data.G′.V′[1:end]] # return set of vertices


# EUC_2D distance
function EUC_dist(u::Vertex, v::Vertex)
    x_sq = (v.pos_x - u.pos_x)^2
    y_sq = (v.pos_y - u.pos_y)^2
    return sqrt(x_sq + y_sq)
    #return floor(sqrt(x_sq + y_sq) + 0.5)
end

contains(p, s) = findnext(s, p, 1) != nothing


function readOPHSData(app::Dict{String,Any})

    str = Unicode.normalize(read(app["instance"], String); stripcc=true)
    breaks_in = [' '; ':'; '\n']
    aux = split(str, breaks_in; limit=0, keepempty=false)
    
    # Initializes data reading
    G′ = InputGraph([], [], Dict())
    data = DataOPHS(G′, [], [], 0.0, 0, [],true,0) # Initializes instance
    data.D = parse(Int, aux[3]) # Number of trips
    data.Tmax = parse(Float64, aux[4]) # Maximum tour length
    
    #= 
    In the input file, aux[2] does not take into account the starting and ending Hotels
    So, the number of hotels is aux[2] + 2 
    The total number of vertices (hotels + customers) = aux[1] + aux[2]
    Set of hotels = {1, 2, ..., aux[2]+2}
    Set of customers = {aux[2]+2+1, ..., aux[1]+aux[2]}
    =#
    data.Hotels = [i for i=1:parse(Int, aux[2])+2]
    data.Customers = [i for i=data.Hotels[end]+1:parse(Int, aux[1])+parse(Int, aux[2])]
    
    pos = 5 # Position initialized at 5 to get remaining info    
    
    # Length of each trip
    for i=1:data.D
        push!(data.Td, parse(Float64, aux[pos])+0.0001)#length of each trip
        pos+=1
    end

    # Is the instance is symmetric? If no, then: data.symmetric = false
    for i=1:length(data.Td)-1,j=i+1:length(data.Td)
        if data.Td[i] != data.Td[j]
            data.symmetric = false
            break
        end
    end

    #data.symmetric = false

    # Getting vertices
    for i=1:parse(Int, aux[1])+parse(Int, aux[2])
        v = Vertex(0, 0, 0, 0) 
        v.id_vertex = i
        v.pos_x = parse(Float64, aux[pos])
        v.pos_y = parse(Float64, aux[pos+1])
        v.score = parse(Float64, aux[pos+2])
        push!(G′.V′, v) # add v in the vertex array
        pos += 3
    end

    #Building the edges set E
    for i=1:length(G′.V′)
        for j=1:length(G′.V′)
            if i < j
                e = (i,j)
                push!(G′.E, e) # add edge e
                data.G′.cost[e] = EUC_dist(data.G′.V′[i],data.G′.V′[j])
            end
        end
    end
    return data
end

edges(data::DataOPHS) = data.G′.E # return set of edges
c(data::DataOPHS,e) = data.G′.cost[e] # cost of the edge e
s(data::DataOPHS,i) = data.G′.V′[i].score # score of the vertex i
d(data::DataOPHS,e) = (e[1] != e[2]) ? data.G′.cost[e] : 0.0 # cost of the edge e
dimension(data::DataOPHS) = length(data.G′.V′) # return number of vertices


# return incident edges of i
function δ(data::DataOPHS, i::Integer)
    incident_edges = Vector{Tuple}()
    for j in 1:i - 1 push!(incident_edges, (j, i)) end
    for j in i + 1:(length(data.G′.V′)) push!(incident_edges, (i, j)) end
    return incident_edges
end
