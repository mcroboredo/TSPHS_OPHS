using VrpSolver, JuMP, ArgParse
using CPLEX

include("data.jl")
include("model.jl")
include("solution.jl")
include("SparseMaxFlowMinCut.jl")



function parse_commandline(args_array::Array{String,1}, appfolder::String)
    s = ArgParseSettings(usage="##### VRPSolver #####\n\n" *
	   "  On interactive mode, call main([\"arg1\", ..., \"argn\"])", exit_after_help=false)
    @add_arg_table s begin
        "instance"
            help = "Instance file path"
        "--cfg", "-c"
            help = "Configuration file path"
            default = "$appfolder/../config/OPHS.cfg"
        "--ub", "-u"
            help = "Upper bound (primal bound)"
            arg_type = Int64
            default = 0
        "--out", "-o"
            help = "Path to write the solution found"
        "--tikz", "-t"
            help = "Path to write the TikZ figure of the solution found."
        "--nosolve", "-n"
            help = "Does not call the VRPSolver. Only to check or draw a given solution."
            action = :store_true
        "--batch", "-b"
            help = "batch file path"
        
        "--complete", "-C"
            help = "is it the complete path generator graph formulation?"
            action = :store_true
    end
   return parse_args(args_array, s)
end

function run_bwtsp(app::Dict{String,Any})
    sol = Solution(0.0, Dict(), [])
    println("Application parameters:")
    for (arg, val) in app
        println("  $arg  =>  $(repr(val))")
    end
    flush(stdout)

    instance_name = split(basename(app["instance"]), ".")[1]

    data = readOPHSData(app)

    solution_found = false
    if !app["nosolve"]
        if data.symmetric
            (model, x, y, g) = build_model(data, app)
        else 
            (model, x, y,t, g) = build_model(data, app)
        end

        optimizer = VrpOptimizer(model, app["cfg"], instance_name)
        set_cutoff!(optimizer, app["ub"]+1.00001)
        
        (status, solution_found) = optimize!(optimizer)

        if solution_found
            if !data.symmetric
                sol = getsolution_non_symmetric(data, optimizer, x, y, t, g, get_objective_value(optimizer), app)
            else 
                sol = getsolution_symmetric(data, optimizer, x, y, g, get_objective_value(optimizer), app)
            end
        end
    end

    println("########################################################")
    if solution_found # Is there a solution?
        #checksolution(data,optimizer,sol, x, y,t,g)
        print_trips(data,sol)
        print_tour(sol)
        println("Cost $(-sol.cost)")
        if data.symmetric
            println("Number of user cuts ", data.cuts)
        end
        
        if app["out"] != nothing
            writesolution(data, app["out"], sol)
        end
        
        if app["tikz"] != nothing
            if data.symmetric
                t = 0
                drawsolution(app["tikz"], data, optimizer, sol,y,t) # write tikz figure
            else 
                drawsolution(app["tikz"], data, optimizer, sol,y,t) # write tikz figure
            end
        end

        #print_stat(data,optimizer,app, g)

    elseif !app["nosolve"]
        if status == :Optimal
            println("Problem infeasible")
        else
            println("Solution not found")
        end
    end

    
    
    println("########################################################")

end

function main(ARGS)
    appfolder = dirname(@__FILE__)
    app = parse_commandline(ARGS, appfolder)
    isnothing(app) && return
    if app["batch"] != nothing
        for line in readlines(app["batch"])
            if isempty(strip(line)) || strip(line)[1] == '#'
                continue
            end
            args_array = [String(s) for s in split(line)]
            app_line = parse_commandline(args_array, appfolder)
            run_bwtsp(app_line)
        end
    else
        run_bwtsp(app)
    end
end

# main()
if isempty(ARGS)
    main(["--help"])
else
    main(ARGS)
end
