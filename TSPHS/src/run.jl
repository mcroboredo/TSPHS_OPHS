using VrpSolver, JuMP, ArgParse

include("data.jl")
include("model.jl")
include("solution.jl")
include("SparseMaxFlowMinCut.jl")

mutable struct statistics
      #name = nome da instância
   #lb_trip = lower bound fornecido pelo usuário para o número de trips
   #opt_trip = último valor fixado para o número de trips.
   #status = string que mostra o status da solução final com a seguinte legenda
      #FIXED -> usuário rodou com um número fixo de trips
      #NOT_FIXED ->  usuário rodou sem um número fixo de trips
      #OPT -> Temos a garantia de que a solução gerada pelo modelo é ótima
      #FEAS -> Temos uma solução mas sem a garantia da otimalidade
      #INF -> Modelo inviável 
      #TL -> Tempo limite foi atingido
   #cutoff = upper bound fornecido pelo usuário. Só está sendo realmente considerado se o número de trips está fixado
   #root_bound = lb no nó raiz do último modelo VRPSolver executado
   #time_bound = tempo no nó raiz do último modelo VRPSolver executado
   #nodes = número de nós do último modelo VRPSolver executado
   #BestDb = Melhor dual bound do último modelo VRPSolver executado
   #BestInc = Melhor solução inteira do último modelo VRPSolver executado
   #last_time = tempo do último modelo VRPSolver executado
   #total_time = tempo total acumulado

   name::String
   lb_trip::Int
   opt_trip::Int
   status::String
   cutoff::Float64
   root_bound::Float64
   root_time::Float64
   nodes::Int
   BestDb::Float64
   BestInc::Float64
   last_time::Float64
   total_time::Float64
   attempts::Array{Int64}
end

function parse_commandline(args_array::Array{String,1}, appfolder::String)
   s = ArgParseSettings(usage="##### VRPSolver arguments #####",exit_after_help=false)
   @add_arg_table s begin
      "instance"
         help = "Instance file path"
      "--cfg", "-c"
         help = "Configuration file path"
         default = "$appfolder/../config/TSPHS.cfg"
      "--ub","-u"
         help = "Upper bound (primal bound)"
         arg_type = Float64
         default = 10000000.0
      "--qvalue","-q"
         help = "minimum number of trips"
         arg_type = Int
         default = 1
      "--Qvalue","-Q"
         help = "maximum number of trips"
         arg_type = Int
         default = 999
      "--sol","-s"
         help = "Solution file path."
      "--out","-o"
         help = "Path to write the solution found"
      "--tikz","-t"
         help = "Path to write the TikZ figure of the solution found."
      "--nosolve","-n"
         help = "Does not call the VRPSolver. Only to check or draw a given solution."
         action = :store_true
      "--batch","-b"
         help = "batch file path"
      "--round","-r"
         help = "We round the distance matrix"
         arg_type = Int
         default = 0
      "--fixedtrip","-f"
         help = "fixed number of trips"
         action = :store_true
      end
   return parse_args(args_array, s)
end

function getconf(time_limit, accumulated_time, maxNodes)
   appfolder = dirname(@__FILE__)
   path_config = string(appfolder, "/../config/TSPHS.cfg")
   path_config_copy = string(appfolder, "/../config/TSPHS_copy.cfg")

   rm(path_config_copy, force = true)
   cp(path_config, path_config_copy)
   #sleep(10)
   
   lines = readlines(path_config_copy)
   
   open(path_config_copy, "w") do f
      for line in lines
         
         if startswith(line, "GlobalTimeLimit")
            new_limit = string(floor(Int, (time_limit - accumulated_time) + 0.5))
            println(f, replace(line, "18000" => new_limit, count = 1))
         else
            println(f, line)
         end
      end
      if maxNodes == true
         println(f, "MaxNbOfBBtreeNodeTreated = 1")
      end
   end

   return path_config_copy
end

function getstatistics(instance_name, q_lb, q, solution_found, accumulated_time, time_limit, rootDB, rootTime,nodes,bestDB,bestINT,last_time,attempts,ub)
   stat = statistics("",0,0,"",0.0,0.0,0.0,0,0.0,0.0,0.0,0.0,[])

   stat.name = instance_name
   stat.lb_trip = q_lb
   stat.opt_trip = q
   if solution_found
      if accumulated_time > time_limit
         stat.status = "NOT_FIXED_FEAS_TL"
      else
         stat.status = "NOT_FIXED_OPT"
      end
   else
      if accumulated_time > time_limit
         stat.status = "NOT_FIXED_INF_TL"
      else
         stat.status = "NOT_FIXED_INF"
      end
   end
   #optimizer.stats[:bcRecRootDb], optimizer.stats[:bcTimeRootEval],optimizer.stats[:bcCountNodeProc],optimizer.stats[:bcRecBestDb],optimizer.stats[:bcRecBestInc]
   #rootDB, rootTime,nodes,bestDB,bestINT
   stat.cutoff = ub
   stat.root_bound = rootDB
   stat.root_time = rootTime
   stat.nodes = nodes
   stat.BestDb = bestDB
   stat.BestInc = bestINT
   stat.last_time = last_time
   stat.total_time = accumulated_time
   stat.attempts = attempts

   return stat
end


function run_TSPHS(app::Dict{String,Any})
   println("Application parameters:")
   for (arg,val) in app
      println("  $arg  =>  $(repr(val))")
   end
   flush(stdout)

   instance_name = split(basename(app["instance"]), ".")[1]
   data = readTSPHSData(app)
   println("########################################################")

   if app["sol"] != nothing
      sol = readsolution(app)
      checksolution(data, sol, app) # checks the solution feasibility
      app["ub"] = (sol.cost < app["ub"]) ? sol.cost : app["ub"] # update the upper bound if necessary
   end

   solution_found = false

   
   
   if !app["nosolve"]

      pos = 0
      for i=1:length(app["instance"]) 
         if app["instance"][i] == '/'
            pos = i
         end
      end
      instance_name = app["instance"][pos + 1:end-4]


      #Get the time limit
      time_limit = 180000 
      str = Unicode.normalize(read(app["cfg"], String); stripcc=true)
      breaks_in = [' '; ':';'='; '\n']
      aux = split(str, breaks_in; limit=0, keepempty=false)
      for i=1:length(aux)
         if aux[i] == "GlobalTimeLimit"
            time_limit = parse(Float64, aux[i+1])
         end
      end
      
      #----------------
      continue_ = true
      accumulated_time = last_time = 0.0
      ub = 10000000.0
      q = app["qvalue"] #minimum number of trips
      time_limit_aux = time_limit
      n_attempts = 6
      Δ = 1.01743287511
      δ = 0.01547247745
      attempts = []
      stat = statistics("",0,0,"",0.0,0.0,0.0,0,0.0,0.0,0.0,0.0,[])
      while continue_
         (model, x) = build_model(data, app, q)
         appfolder = dirname(@__FILE__)
         
         #--------- Root node --------------
         maxNodes = true
         path_config_copy = getconf(time_limit,accumulated_time,maxNodes)
         optimizer = VrpOptimizer(model, path_config_copy, instance_name)
         #if app["fixedtrip"]
            #The initial upper bound only can be used if the number of trips is fixed
            #set_cutoff!(optimizer, app["ub"]+0.1)
         #end
         (status, solution_found) = optimize!(optimizer)
         rm(path_config_copy, force = true)
         accumulated_time += optimizer.stats[:bcTimeMain] / 100
         time_limit_aux -= optimizer.stats[:bcTimeMain] / 100
         if solution_found || accumulated_time >= time_limit
            if solution_found
               println("The solution was found during the root node")
               sol = getsolution(data, optimizer, x, get_objective_value(optimizer), app)
            else
               println("The TL was reached after root node")
            end
            break
         end
         println("end of the root node ", q, " ", optimizer.stats[:bcRecRootDb])
         #sleep(10)
         lb = optimizer.stats[:bcRecRootDb]
         #ub = (Δ+δ)*optimizer.stats[:bcRecRootDb]
         #----------- End root node -------
         #=
         if app["fixedtrip"]
            last_time = optimizer.stats[:bcTimeMain] / 100

            stat.name = instance_name
            stat.lb_trip = app["qvalue"]
            stat.opt_trip = app["qvalue"]
            if solution_found
               if last_time > time_limit
                  stat.status = "FIXED_FEAS_TL"
               else
                  stat.status = "FIXED_OPT"
               end
            else
               if last_time > time_limit
                  stat.status = "FIXED_INF_TL"
               else
                  stat.status = "FIXED_INF"
               end
            end
            stat.cutoff = app["ub"]
            stat.root_bound = optimizer.stats[:bcRecRootDb]
            stat.root_time = optimizer.stats[:bcTimeRootEval]
            stat.nodes = optimizer.stats[:bcCountNodeProc]
            stat.BestDb = optimizer.stats[:bcRecBestDb]
            stat.BestInc = optimizer.stats[:bcRecBestInc]
            stat.last_time = last_time
            stat.total_time = accumulated_time
            break
         end
         =#

         #------------------------
         if optimizer.stats[:bcRecRootDb] < 10000000.0
            for i=1:n_attempts
               if accumulated_time < time_limit
                  maxNodes = false
                  path_config_copy = getconf(time_limit,accumulated_time,maxNodes)
                  optimizer = VrpOptimizer(model, path_config_copy, instance_name)
                  ub = lb*(Δ+i*δ)
                  set_cutoff!(optimizer, ub)
                  (status, solution_found) = optimize!(optimizer)
                  rm(path_config_copy, force = true)
                  accumulated_time += optimizer.stats[:bcTimeMain] / 100
                  time_limit_aux -= optimizer.stats[:bcTimeMain] / 100
                  if solution_found || accumulated_time >= time_limit
                     push!(attempts,i)
                     last_time = optimizer.stats[:bcTimeMain] / 100
                     
                     stat = getstatistics(instance_name, app["qvalue"], q, solution_found, accumulated_time, time_limit, optimizer.stats[:bcRecRootDb], optimizer.stats[:bcTimeRootEval],optimizer.stats[:bcCountNodeProc],optimizer.stats[:bcRecBestDb],optimizer.stats[:bcRecBestInc],last_time,attempts,ub)
                     
                     if solution_found
                        println("The solution was found during the attempt ", i)
                        sol = getsolution(data, optimizer, x, get_objective_value(optimizer), app)
                     else
                        println("The TL was reached after attempt ", i)
                     end
                     continue_ = false
                     break
                  end
                  println("end of the attempt ", i, " ", q, " ", ub)
                  #sleep(10)
                  
               end
            end
            #---------------------------

            if continue_ == false
               break
            end

            #----- Running without upper bound -----
            maxNodes = false
            path_config_copy = getconf(time_limit,accumulated_time,maxNodes)
            optimizer = VrpOptimizer(model, path_config_copy, instance_name)
            #set_cutoff!(optimizer, ub)
            (status, solution_found) = optimize!(optimizer)
            rm(path_config_copy, force = true)
            accumulated_time += optimizer.stats[:bcTimeMain] / 100
            time_limit_aux -= optimizer.stats[:bcTimeMain] / 100

            if solution_found || accumulated_time >= time_limit
               ub = 10000000000.0
               stat = getstatistics(instance_name, app["qvalue"], q, solution_found, accumulated_time, time_limit, optimizer.stats[:bcRecRootDb], optimizer.stats[:bcTimeRootEval],optimizer.stats[:bcCountNodeProc],optimizer.stats[:bcRecBestDb],optimizer.stats[:bcRecBestInc],last_time,attempts, ub)
               if solution_found
                  println("The solution was found without upper bound")
                  sol = getsolution(data, optimizer, x, get_objective_value(optimizer), app)
               else
                  println("The TL was reached after running without upper bound")
               end
               break
            end

            println("end of running without upper bound ", q)
            push!(attempts,n_attempts + 1)
            #sleep(10)
            #-------------
         
         else
            println("the problem was proved to be infeasible for ", q, " trips during the root node")
            push!(attempts,0)
         
         
         end

         
         q = q + 1
         println("Nem value of q ", q)
         #sleep(5)

         #=
         if feasible == false && accumulated_time <= time_limit
            q += 1
         else
            last_time = optimizer.stats[:bcTimeMain] / 100

            stat.name = instance_name
            stat.lb_trip = app["qvalue"]
            stat.opt_trip = q
            if solution_found
               if last_time > time_limit
                  stat.status = "NOT_FIXED_FEAS_TL"
               else
                  stat.status = "NOT_FIXED_OPT"
               end
            else
               if last_time > time_limit
                  stat.status = "NOT_FIXED_INF_TL"
               else
                  stat.status = "NOT_FIXED_INF"
               end
            end
            stat.cutoff = 10000000.0
            stat.root_bound = optimizer.stats[:bcRecRootDb]
            stat.root_time = optimizer.stats[:bcTimeRootEval]
            stat.nodes = optimizer.stats[:bcCountNodeProc]
            stat.BestDb = optimizer.stats[:bcRecBestDb]
            stat.BestInc = optimizer.stats[:bcRecBestInc]
            stat.last_time = last_time
            stat.total_time = accumulated_time
         end
         =#

      end
   
   end
   println(" ")

   println("######### Customized statistics #########################")

   println("Customized_statistics: ", stat.name, " ", stat.lb_trip, " ",stat.opt_trip, " ",stat.status, " ",stat.cutoff, " ",stat.root_bound, " ",stat.root_time, " ",stat.nodes, " ",stat.BestDb, " ",stat.BestInc, " ",stat.last_time, " ",stat.total_time, " ", stat.attempts)
   
   println("########################################################")
   println("")
   println("########################################################")
   
   if solution_found || app["sol"] != nothing # Is there a solution?
      checksolution(data, sol, app)
      print_routes(data, sol)
      println("Cost $(sol.cost)")
      if app["out"] != nothing
         writesolution(app["out"], sol)
      end
      if app["tikz"] != nothing
         drawsolution(app["tikz"], data, sol) # write tikz figure
      end
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
         run_TSPHS(app_line)
      end
   else
      run_TSPHS(app)
   end
end

# main()
if isempty(ARGS)
   main(["--help"])
else
   main(ARGS)
end
