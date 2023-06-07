import Unicode

function distance(arc::Tuple{Int64, Int64}, coord_x,coord_y,st)
    x_sq = (coord_x[arc[1]] - coord_x[arc[2]])^2
    y_sq = (coord_y[arc[1]] - coord_y[arc[2]])^2
    return floor(sqrt(x_sq + y_sq) + 0.5)
    return sqrt(x_sq + y_sq) + 0.5*(st[arc[1]]+st[arc[2]])
end

include("concorde_tsp.jl")


function main()
    str = Unicode.normalize(read("data/berlin52.tsp", String); stripcc=true)
    breaks_in = [' '; ':'; '\n']
    aux = split(str, breaks_in; limit=0, keepempty=false)

    n = parse(Int, aux[1])
    @show n
    coord_x = [  ]
    coord_y = [  ]
    st = [  ]

    pos = 2
    for i=1:n
        push!(coord_x,parse(Float64,aux[pos]))
        push!(coord_y,parse(Float64,aux[pos+1]))
        push!(st,parse(Float64,aux[pos+2]))
        pos += 3
    end

    
    @show coord_x
    @show coord_y
    @show st

    cost = zeros(Float64,n,n)
    for i=1:n
        for j=1:n
            cost[i,j] = distance((i,j),coord_x,coord_y,st)
        end
    end


    matriz = rand(5,5)
    cost_, tour = concorde_st(5,matriz)
    @show cost_, tour
   

   

end

main()


