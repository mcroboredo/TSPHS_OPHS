include("extra/factorials.jl")
include("extra/permutations.jl")


mutable struct Solution
    cost::Union{Int,Float64}
    trips::Dict{Float64,Array{Int64}}
    tour::Array{Int}
end

mutable struct Custom_Solution
    ub0::Float64
    opt::Float64
    root_lb::Float64
    nodes::Int64
    cuts::Int64
    time::Float64
    avg_customers::Float64
end

function convert_time(time)
    if time < 100
        return 0.0
    end 
    a = string(time)
    b = string(a[1:length(a)-2],".",a[length(a)-1:length(a)])
    return parse(Float64, b)
end

# build Solution from the variables x
function getsolution_non_symmetric(data::DataOPHS, optimizer::VrpOptimizer, x, y,t,g,objval, app::Dict{String,Any})
    
    D = data.D # Number of trips
    T =[i for i=1:D]
    
    sol_trips, sol_edges, sol_cost = Dict(),[], get_objective_value(optimizer)

    #getting each trip
    for path_id in 1:get_number_of_positive_paths(optimizer)
        trip, trip_edges, source, sink, id = [], [], 0,0,0
        
        #getting the hotel where the trip begins, the hotel where the trip ends and the trip index
        for k in T
            for i in data.Hotels
                val2 = get_value(optimizer, y[i,k], path_id)
                if val2 > 0.99
                    id = k
                    source = i
                end

                val3 = get_value(optimizer, t[i,k], path_id)
                if val3 > 0.99
                    id = k
                    sink = i
                end
            end
        end
        #getting the edges of each trip
        for e in edges(data)
            val = get_value(optimizer, x[e], path_id)
            if val > 0
                push!(trip_edges, e)
                push!(sol_edges, e)
            end
        end
        #getting the complete trip
        push!(trip, source)
        next = source
        visited = Dict()
        for e in trip_edges
            visited[e] = false
        end
        for e in trip_edges
            if e[1] == source || e[2] == source
                visited[e] = true
                if e[1] == source
                    next = e[2]
                    push!(trip, next)
                else
                    next = e[1]
                    push!(trip, next)
                end
                break
            end
        end
        finish = false
        while finish == false
            finish = true
            for e in trip_edges
                if visited[e] == false
                    if e[1] == next || e[2] == next
                        finish = false
                        visited[e] = true
                        if e[1] == next
                            next = e[2]
                            push!(trip, next)
                        else
                            next = e[1]
                            push!(trip, next)
                        end
                    end
                end
            end
        end
        sol_trips[id] = trip
    end

    tour = []
    for k=1:D
        for j=1:length(sol_trips[k])-1
            push!(tour,sol_trips[k][j])
        end
    end
    push!(tour, 2)
    return Solution(sol_cost, sol_trips,tour)
end

function getsolution_symmetric(data::DataOPHS, optimizer::VrpOptimizer, x, y,g,objval, app::Dict{String,Any})
    id = 1
    D = data.D # Number of trips
    T =[i for i=1:D]
    
    sol_trips, sol_edges, sol_cost = Dict(),[], get_objective_value(optimizer)

    #getting each trip
    for path_id in 1:get_number_of_positive_paths(optimizer)
        trip, trip_edges, source, sink= [], [], 0,0
        #getting the 2 hotels in the trip
        hs = 0
        for i in data.Hotels
            val = get_value(optimizer, y[i], path_id)
            if val > 0
                hs = i 
            end
        end
        #getting the edges of each trip
        for e in edges(data)
            val = get_value(optimizer, x[e], path_id)
            if val > 0
                push!(trip_edges, e)
            end
        end
        if !((1,2) in trip_edges)
             #getting the complete trip
            push!(trip, hs)
            next = hs
            visited = Dict()
            for e in trip_edges
                visited[e] = false
            end
            for e in trip_edges
                if e[1] == hs || e[2] == hs
                    visited[e] = true
                    if e[1] == hs
                        next = e[2]
                        push!(trip, next)
                    else
                        next = e[1]
                        push!(trip, next)
                    end
                    break
                end
            end
            finish = false
            while finish == false
                finish = true
                for e in trip_edges
                    if visited[e] == false
                        if e[1] == next || e[2] == next
                            finish = false
                            visited[e] = true
                            if e[1] == next
                                next = e[2]
                                push!(trip, next)
                            else
                                next = e[1]
                                push!(trip, next)
                            end
                        end
                    end
                end
            end
            sol_trips[id] = trip
            id += 1
        end
    end
    tour = []
    trips_array = []
    for i=1:D
        push!(trips_array, sol_trips[i])
    end

    T = collect(permutations(trips_array))
    found_tour = false
    for i=1:length(T)
        found_tour = true
        if T[i][1][1] == 1 && T[i][D][end] == 2
            for k=1:D-1
                if T[i][k][end] != T[i][k+1][1]
                    found_tour = false
                end
            end
        else
            found_tour = false
        end

        if found_tour
            for k=1:D
                for j=1:length(T[i][k])-1
                    push!(tour,T[i][k][j])
                end
            end
            push!(tour, 2)
            @show tour
            break
        end
    end

    


    return Solution(sol_cost, sol_trips,tour)
end

function print_trips(data::DataOPHS, solution)
    for k=1:length(solution.trips)
        print("trip ", k, " : ")
        for i=1:length(solution.trips[k])
            print(solution.trips[k][i], " ")
        end
        println()
    end
end

function print_tour(solution)
    print("tour : ")
    for k=1:length(solution.tour)
        if k == length(solution.tour)
            print(solution.tour[k])
        else
            print(solution.tour[k], " ")
        end
    end
    println()
end

# checks the feasiblity of a solution
ed(i, j) = i < j ? (i, j) : (j, i)
function checksolution(data::DataOPHS, optimizer::VrpOptimizer, solution, x, y,t,g)

    if data.symmetric
        #checking the number of trips of the solution 
        if get_number_of_positive_paths(optimizer) != data.D + 1
            println("error: number of trips of the solution considering the dummy trip = ", get_number_of_positive_paths(optimizer), " Instance number of trips  considering the dummy trip = ", data.D+1)
        end
    else
        #checking the number of trips of the solution 
        if get_number_of_positive_paths(optimizer) != data.D
            println("error: number of trips of the solution = ", get_number_of_positive_paths(optimizer), " Instance number of trips  = ", data.D)
        end
    end

    #cheking the objective function
    objective = 0.0
    for c in data.Customers
        val = get_value(optimizer, g[c])
        if val > 0.99
            objective += s(data,c)
        end
    end
    if -objective < get_objective_value(optimizer) - 0.001 || -objective > get_objective_value(optimizer) + 0.001
        println("error: objective calculated by hand ", -objective, " objective calculated by VRPSolver ", get_objective_value(optimizer))
    end

    
    if !data.symmetric
        #checking if the tour starts at the hotel 1
        if get_value(optimizer, y[1,1]) != 1
            println("error: the tour does not start at the hotel 1")
        end

        #checking if the tour ends at the hotel 2
        if get_value(optimizer, t[2,data.D]) != 1
            println("error: the tour does not end at the hotel 2")
        end

        #cheking if y[i,k] = t[i,k-1]
        for k=2:data.D, i in data.Hotels
            if get_value(optimizer, y[i,k]) != get_value(optimizer, t[i,k-1])
                println("error: y[",i,",",k,"] = ", get_value(optimizer, y[i,k]), " t[",i,",",k-1,"] = ", get_value(optimizer, t[i,k-1]))
            end
        end

        #cheking if each trip starts at one hotel exactly
        for k=1:data.D
            sum = 0
            for i in data.Hotels
                sum+=get_value(optimizer, y[i,k])
            end
            if sum != 1
                println("error: the trip ", k, " starts at ", sum, " hotels")
            end
        end

        #cheking if each trip ends at one hotel exactly
        for k=1:data.D
            sum = 0
            for i in data.Hotels
                sum+=get_value(optimizer, t[i,k])
            end
            if sum != 1
                println("error: the trip ", k, " ends at ", sum, " hotels")
            end
        end
    end

    #cheking if the maximum tour length is being respected
    total_tour_length = 0.0
    for e in solution.edges
        total_tour_length += d(data,e)
    end 
    if total_tour_length > data.Tmax
        println("error: the maximum tour length is not being respected: ")
        println("length of the solution = ", total_tour_length, " total length of the instance = ", data.Tmax)
    end

    #checking if each maximum trip length is being respected 
    for k=1:data.D
        trip, trip_len = solution.trips[k], 0.0
        for i=1:length(trip)-1
            #e = (trip[i], trip[i+1])
            trip_len+=d(data,ed(trip[i], trip[i+1]))
        end
        if trip_len > data.Td[k]
            println("error: the maximum length of the trip ", k, " is not being respected: ")
            println("length of the trip = ", trip_len, " total trip length of the instance = ", data.Td[k])
        end
    end
end

# write solution in a file
function writesolution(data::DataOPHS, solpath, solution)
    open(solpath, "w") do f
        for k=1:data.D
            write(f, "trip $(k) :")
            for i=1:length(solution.trips[k])
                write(f, " $(solution.trips[k][i])")
            end
            write(f, "\n")
        end
        for k=1:length(solution.tour)
            write(f, "tour :")
            if k == length(solution.tour)
                write(f,"$(solution.tour[k])\n")
            else
                write(f,"$(solution.tour[k]) ")
            end
        end
        write(f, "Cost : $(-solution.cost)\n")
    end
end

# write solution as TikZ figure (.tex)
function drawsolution(tikzpath,data::DataOPHS, optimizer::VrpOptimizer, solution,y,t)
    open(tikzpath, "w") do f
        write(f, "\\documentclass[crop,tikz]{standalone}\n\\begin{document}\n")
        # get limits to draw
        pos_x_vals = [i.pos_x for i in data.G′.V′]
        pos_y_vals = [i.pos_y for i in data.G′.V′]

        #if data.insType == "EUC_2D" || data.insType == "GEO"
        scale_fac = 1 / (max(maximum(pos_x_vals), maximum(pos_y_vals)) / 10)
        write(f, "\\begin{tikzpicture}[thick, scale=1, every node/.style={scale=0.3}]\n")
        for i in data.G′.V′
            if i.id_vertex in data.Customers
                x_plot = scale_fac * i.pos_x
                y_plot = scale_fac * i.pos_y
                write(f, "\t\\node[draw, line width=0.1mm, circle, fill=white, inner sep=0.05cm] (v$(i.id_vertex)) at ($(x_plot),$(y_plot)) {\\footnotesize $(i.id_vertex)};\n")
            end
        end

        for i in data.G′.V′
            if i.id_vertex in data.Hotels
                x_plot = scale_fac * i.pos_x
                y_plot = scale_fac * i.pos_y
                sum = 0
                if !data.symmetric
                    for k=1:data.D
                        sum += get_value(optimizer, y[i.id_vertex,k]) + get_value(optimizer, t[i.id_vertex,k])
                    end
                else
                    sum += get_value(optimizer, y[i.id_vertex])
                end
                if sum < 1 #the hotel is not used
                    write(f, "\t\\node[draw, line width=0.1mm, circle, fill=black, inner sep=0.05cm, text=white] (v$(i.id_vertex)) at ($(x_plot),$(y_plot)) {\\footnotesize $(i.id_vertex)};\n")
                else
                    write(f, "\t\\node[draw, line width=0.1mm, circle, fill=red, inner sep=0.05cm, text=white] (v$(i.id_vertex)) at ($(x_plot),$(y_plot)) {\\footnotesize $(i.id_vertex)};\n")
                end
            end
        end
        edge_style = "-,line width=0.8pt"
        for k=1:data.D
            for i=1:length(solution.trips[k])-1
                e_1,e_2 = solution.trips[k][i], solution.trips[k][i+1]
                if solution.trips[k][i] > solution.trips[k][i+1]
                    e_1,e_2 = solution.trips[k][i+1], solution.trips[k][i]
                end
                write(f, "\t\\draw[$(edge_style)] (v$(e_1)) -- (v$(e_2));\n")
            end
        end
        write(f, "\\end{tikzpicture}\n")
        write(f, "\\end{document}\n")
    end
end

function print_stat(data::DataOPHS, optimizer::VrpOptimizer,app::Dict{String,Any}, g)
    instance_name = split(basename(app["instance"]), ".")[1]
    
    println("########################################################")

    println("columns: name |C| |H| Tmax :bcRecRootDb :bcTimeRootEval :bcCountNodeProc :bcRecBestDb :bcTimeMain :bcCountMastSol :bcCountCol :bcCountCutInMaster :bcTimeMastMPsol :bcTimeCgSpOracle :bcTimeCutSeparation :bcTimeAddCutToMaster :bcTimeSetMast :bcTimeRedCostFixAndEnum :bcTimeEnumMPsol :bcTimeSBphase1 :bcTimeSBphase2 :bcTimePrimalHeur :bcRecBestInc")
    println("complete_stat: ", instance_name, " ", length(data.Customers), " ", length(data.Hotels), " ", data.Tmax, " ", optimizer.stats[:bcRecRootDb], " ", convert_time(optimizer.stats[:bcTimeRootEval]), " ", optimizer.stats[:bcCountNodeProc], " ",
    optimizer.stats[:bcRecBestDb], " ", convert_time(optimizer.stats[:bcTimeMain]), " ", optimizer.stats[:bcCountMastSol], " ",
    optimizer.stats[:bcCountCol], " ", optimizer.stats[:bcCountCutInMaster], " ", convert_time(optimizer.stats[:bcTimeMastMPsol]), " ",
    convert_time(optimizer.stats[:bcTimeCgSpOracle]), " ", convert_time(optimizer.stats[:bcTimeCutSeparation]), " ", convert_time(optimizer.stats[:bcTimeAddCutToMaster]), " ",
    convert_time(optimizer.stats[:bcTimeSetMast]), " ", convert_time(optimizer.stats[:bcTimeRedCostFixAndEnum]), " ", convert_time(optimizer.stats[:bcTimeEnumMPsol]), " ",
    convert_time(optimizer.stats[:bcTimeSBphase1]), " ", convert_time(optimizer.stats[:bcTimeSBphase2]), " ", convert_time(optimizer.stats[:bcTimePrimalHeur]), " ",optimizer.stats[:bcRecBestInc])

    println("########################################################")

    number_of_customers = 0
    for c in data.Customers
        if get_value(optimizer, g[c]) > 0.0
            number_of_customers += 1
        end
    end
    average_customers_trip = number_of_customers/(get_number_of_positive_paths(optimizer))
    time = convert_time(optimizer.stats[:bcTimeMain])
    sol = Custom_Solution(app["ub"],optimizer.stats[:bcRecBestInc], optimizer.stats[:bcRecRootDb], trunc(optimizer.stats[:bcCountNodeProc]), data.cuts, time, average_customers_trip )
    
    println("cols name |C| |H| Tmax ub0 opt avg_customers_by_route lb0 nodes user_cuts time")
    println("custom_stat: ", instance_name, " ", length(data.Customers), " ", length(data.Hotels), " ", data.Tmax, " ", sol.ub0, " ", sol.opt, " ", sol.avg_customers, " ", sol.root_lb, " ", sol.nodes, " ", sol.cuts, " ", sol.time)

    println("########################################################")
end

