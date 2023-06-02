mutable struct Solution
   cost::Union{Int,Float64}
   routes::Array{Array{Int}}
   edge_sol::Dict{Tuple{Int64,Int64},Float64}
end



# build Solution from the variables x
function getsolution_alt(data::DataTSPHS, optimizer::VrpOptimizer, x, objval, app::Dict{String,Any})
   E, dim, H = edges(data), dimension(data), data.H′
   adj_list = [[] for i in 1:dim]
   incidence = [0 for i in H]
   edge_solution = Dict()

   trips = []
   cont = 1
   for path_id in 1:get_number_of_positive_paths(optimizer)
      trip, trip_edges, source, sink= [], [], 0,0

      println("Trip ", cont)
   
      #getting the edges of each trip
      for e in E
          val = get_value(optimizer, x[e], path_id)
          if val > 0
              @show e
              push!(trip_edges, e)
          end
      end

      current = 0
      visited = Dict()
      for e in trip_edges
         visited[e] = false
      end
      for e in trip_edges
         if e[1] in H
            @show e
            current = e[2]
            visited[e] = true
            trip=[e[1],e[2]]
            break
         end
      end

      continue_ = true

      while continue_
         for e in trip_edges
            if visited[e] == false && current in e
               @show e
               visited[e] = true
               if e[1] == current
                  current = e[2]
                  push!(trip,current)
               else
                  current = e[1]
                  push!(trip,current)
               end
               break
            end
         end
         continue_ = false
         for e in trip_edges
            if visited[e] == false
               continue_ = true
               break
            end
         end

         
      end
      push!(trips, trip)

      cont+= 1
      @show trip
   end

   
   #return Solution(objval, r_aux, edge_solution)
end

# build Solution from the variables x
function getsolution(data::DataTSPHS, optimizer::VrpOptimizer, x, objval, app::Dict{String,Any})
   E, dim, H = edges(data), dimension(data), data.H′
   adj_list = [[] for i in 1:dim]
   incidence = [0 for i in H]
   edge_solution = Dict()

   for (i,j) in E
      e = (i,j)
      val = get_value(optimizer, x[e])
      if val > 0.5
         push!(adj_list[e[1]], e[2])
         push!(adj_list[e[2]], e[1])
         if in(e[1],H) incidence[e[1]]+=1 end
         if in(e[2],H) incidence[e[2]]+=1 end
         if val > 1.5
            push!(adj_list[e[1]], e[2])
            push!(adj_list[e[2]], e[1])
            if in(e[1],H) incidence[e[1]]+=1 end
            if in(e[2],H) incidence[e[2]]+=1 end
         end
         edge_solution[e] = val
      end
   end
   r_aux = []
   visited, routes = [false for i in 1:dim], [[] for i in H]
   tour_check = [0 for i in H]
   for h in 1:length(H)
      if !isempty(adj_list[h])
         h_aux = h
         for i in adj_list[h_aux]
            r = []
            if in(i,H)
               push!(r,h_aux)
               push!(r,i)
               visited[i] = true
            elseif !visited[i]
               prev = h_aux
               push!(r, h_aux)
               push!(r, i)
               visited[i] = true
               filter!(a->a!=prev,adj_list[i])
               next, prev = !isempty(adj_list[i]) ? adj_list[i][1] : prev , i
               while !in(next,H)
                  push!(r, next)
                  visited[next] = true
                  aux = next
                  filter!(a->a!=prev,adj_list[next])
                  next, prev = !isempty(adj_list[next]) ? adj_list[next][1] : prev , aux
               end
               push!(r, next)
               visited[next] = true
               r[end] == h_aux ? pushfirst!(routes[h_aux],r) : push!(routes[h_aux],r)
               r[end] == h_aux ? tour_check[h_aux] += 1 : nothing
               push!(r_aux,r)
            end
         end
      else
         visited[h] = true
      end
   end

   H′ = []
   for i in 1:length(incidence)
      incidence[i] > 2 ? push!(H′, (i, incidence[i])) : nothing
   end
   sort!(H′, by = x -> x[2])
   # println(H′)
   # println(incidence)
   for i in H
      if incidence[i] == 2
         r1, del = [], []
         for j in 1:length(r_aux)
            flag = 0
            if r_aux[j][1] == i
               isempty(r1) ? reverse!(r_aux[j]) : popfirst!(r_aux[j])
               append!(r1,r_aux[j])
               push!(del,r_aux[j])
               flag+=1
            elseif r_aux[j][end] == i
               isempty(r1) ? nothing : popfirst!(reverse!(r_aux[j]))
               append!(r1,r_aux[j])
               push!(del,r_aux[j])
               flag+=1
            end
            flag == 2 ? break : nothing
         end
         setdiff!(r_aux,del)
         r1[1] == r1[end] ? tour_check[r1[1]]+=1 : nothing
         r1[1] == r1[end] ? pushfirst!(r_aux,r1) : push!(r_aux,r1)
         incidence[i] = 0
      end
   end
   # println(r_aux)
   # println()
   for k in H′
      i = k[1]
      if tour_check[i] > 0
         r1, del, aux = [], [], []
         for j in 1:length(r_aux)
            flag = 0
            if r_aux[j][1] == i && r_aux[j][end] == i
               isempty(r1) ? nothing : popfirst!(r_aux[j])
               append!(r1,r_aux[j])
               push!(del,r_aux[j])
               flag+=1
            elseif (r_aux[j][1] == i || r_aux[j][end] == i) && length(aux) < 2
               push!(aux,r_aux[j])
               push!(del,r_aux[j])
            end
            (flag == tour_check[i] && length(aux) == 2) ? break : nothing
         end
         if !isempty(aux)
            aux2 = []
            aux[1][1] == i ? pop!(reverse!(aux[1])) : pop!(aux[1])
            aux[2][1] == i ? popfirst!(aux[2]) : popfirst!(reverse!(aux[2]))
            append!(aux2,aux[1])
            append!(aux2,r1)
            append!(aux2,aux[2])
            aux2[1] == aux2[end] ? pushfirst!(r_aux,aux2) : push!(r_aux,aux2)
            aux2[1] == aux2[end] ? tour_check[aux2[1]]+=1 : nothing
         else
            r1[1] == r1[end] ? pushfirst!(r_aux,r1) : push!(r_aux,r1)
            r1[1] == r1[end] ? tour_check[r1[1]]+=1 : nothing
         end
         setdiff!(r_aux,del)
         incidence[i], tour_check[i] = incidence[i]-2*tour_check[i]-length(aux), 0
      end
      length(r_aux) <= 2 ? break : nothing
   end
   # println(r_aux)
   # println()
   i = 2
   while length(r_aux) != 1
      flag = false
      if r_aux[i][1] == r_aux[1][end]
         popfirst!(r_aux[i])
         flag = true
      elseif r_aux[i][end] == r_aux[1][end]
         popfirst!(reverse!(r_aux[i]))
         flag = true
      end
      if flag
         r1, del = [], []
         append!(r1,r_aux[1])
         push!(del,r_aux[1])
         append!(r1,r_aux[i])
         push!(del,r_aux[i])
         r1[1] == r1[end] ? pushfirst!(r_aux,r1) : push!(r_aux,r1)
         setdiff!(r_aux,del)
         i = 1
      end
      i+=1
   end

   objval = round(objval, digits=2)
   return Solution(objval, r_aux, edge_solution)
end

function print_routes(data::DataTSPHS, solution)
   for (i,r) in enumerate(solution.routes)
      print("Tour #$i: ")
      for j in r
         !in(j,data.H′) ? print("$(j-1) ") : printstyled("$(j-1) ", color=:cyan)
      end
      println()
   end
end

# checks the feasiblity of a solution
function checksolution(data::DataTSPHS, solution, app)
   dim, L = dimension(data), data.Lim
   visits = [0 for i in 1:dim]
   sum_cost = 0.0
   for (i,r) in enumerate(solution.routes)
      sum_time, prev = 0.0, r[1]
      for j in r[2:end]
         visits[j] += 1
         (visits[j] == 2 && in(j,data.C′)) && error("Customer $j was visited more than once")
         sum_cost += distance(data, (prev,j))
         in(j,data.C′) ? sum_time += (service(data,j) + distance(data, (prev,j))) : sum_time += distance(data, (prev,j))
         (sum_time > L +0.001) && error("Route is violating the daily time limit. The daily time limit is $(sum_time) and L is $(L)")
         !in(j,data.C′) ? sum_time = 0 : nothing
         prev = j
      end
   end
   !isempty(filter(a->a==0,visits[nb_hotels(data)+1:end])) && error("The following customers were not visited: $(filter(a->a==0,visits[nb_hotels(data)+1:end]))")
   (abs(solution.cost-sum_cost) > 0.001) && error("Cost calculated from the routes ($sum_cost) is different from that passed as"*
                                                                                        " argument ($(solution.cost)).")
end
# read solution from file (CVRPLIB format)
contains(p, s) = findnext(s, p, 1) != nothing
function readsolution(app::Dict{String,Any})
   str = read(app["sol"], String)
   breaks_in = [' '; ':'; '\n';'\t';'\r']
   aux = split(str, breaks_in; limit=0, keepempty=false)
   sol = Solution(0, [], Dict())
   j = 3
   while j <= length(aux)
      r = []
      while j <= length(aux)
         push!(r, parse(Int, aux[j])+1)
         j += 1
         if contains(lowercase(aux[j]), "cost") || contains(lowercase(aux[j]), "route")
            break
         end
      end
      push!(sol.routes, r)
      if contains(lowercase(aux[j]), "cost")
         sol.cost = parse(Float64, aux[j+1])
         return sol
      end
      j += 2 # skip "Route" and "#j:" elements
   end
   error("The solution file was not read successfully.")
   return sol
end

# write solution in a file
function writesolution(solpath, solution)
   open(solpath, "w") do f
      for (i,r) in enumerate(solution.routes)
         write(f, "Route #$i: ")
         for j in r
            write(f, "$(j-1) ")
         end
         write(f, "\n")
      end
      write(f, "Cost $(solution.cost)\n")
   end
end

# write solution as TikZ figure (.tex)
function drawsolution(tikzpath, data, solution)
   open(tikzpath, "w") do f
      write(f,"\\documentclass[crop,tikz]{standalone}\n\\usetikzlibrary{shapes.geometric}\n\\begin{document}\n")
      # get limits to draw
      pos_x_vals = [i.pos_x for i in data.G′.V′]
      pos_y_vals = [i.pos_y for i in data.G′.V′]

      scale_fac = 1/(max(maximum(pos_x_vals),maximum(pos_y_vals))/10)

      write(f,"\\begin{tikzpicture}[thick, scale=1, every node/.style={scale=0.3},square/.style={regular polygon,regular polygon sides=4}]\n")
      for i in data.G′.V′
         x_plot = scale_fac*i.pos_x
         y_plot = scale_fac*i.pos_y
         if i.id_vertex in data.H′ # plot balcks
            fill = ""
            i.id_vertex == data.H′[1] ? fill = "cyan" : fill = "yellow"
            write(f, "\t\\node[draw, line width=0.1mm, square, fill=$fill, inner sep=0.05cm] (v$(i.id_vertex)) at ($(x_plot),$(y_plot)) {\\footnotesize $(i.id_vertex-1)};\n")
         else
            write(f, "\t\\node[draw, line width=0.1mm, circle, fill=white, inner sep=0.05cm] (v$(i.id_vertex)) at ($(x_plot),$(y_plot)) {\\footnotesize $(i.id_vertex-1)};\n")
         end
      end

      if !isempty(solution.edge_sol)
         for e in keys(solution.edge_sol)
            edge_style = "-,line width=0.2mm"
            if solution.edge_sol[e] > 1.5
               write(f, "\t\\draw[$(edge_style),bend right=15] (v$(e[1])) to (v$(e[2]));\n")
               write(f, "\t\\draw[$(edge_style),bend left=15] (v$(e[1])) to (v$(e[2]));\n")
            else
               write(f, "\t\\draw[$(edge_style)] (v$(e[1])) to (v$(e[2]));\n")
            end
         end
      else
         aux = []
         for r in solution.routes
            prev = r[1]
            edge_style = "-,line width=0.2mm"
            for i in r[2:end]
               e = prev < i ? (prev,i) : (i,prev)
               if in(e,aux)
                  write(f, "\t\\draw[$(edge_style),bend right=15] (v$(e[1])) to (v$(e[2]));\n")
               else
                  write(f, "\t\\draw[$(edge_style)] (v$(e[1])) to (v$(e[2]));\n")
                  push!(aux, e)
               end
               prev = i
            end
         end
      end
      write(f, "\\end{tikzpicture}\n")
      write(f, "\\end{document}\n")
   end
end
