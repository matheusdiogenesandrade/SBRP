module Main

include("data.jl")
include("model.jl")
include("solve.jl")
include("sol.jl")

using .Data
using .Data.SBRP
using .Model
using .Model.ModelSBRP
using .Model.ModelSBRPComplete
using .Model.ModelTSP
using .Model.ModelSBRPMax
using .Solution
using ArgParse

function parse_commandline(args_array::Array{String,1}, appfolder::String)
  s = ArgParseSettings(usage="  On interactive mode, call main([\"arg1\", ..., \"argn\"])", exit_after_help=false)
  @add_arg_table s begin
    "instance"
    help = "Instance file path"
    "--instance-type"
    help = "Instance type (carlos, matheus)"
    default = "carlos"
    "--nosolve","-n"
    help = "Not solve flag"
    action = :store_true
    "--out","-o"
    help = "Path to write the solution found"
    "--batch","-b" 
    help = "batch file path" 
  end
  return parse_args(args_array, s)
end
 
function run(app::Dict{String,Any})
  println("Application parameters:")
  for (arg,val) in app
    println("  $arg  =>  $(repr(val))")
  end
  flush(stdout)
  # read instance
  readInstanceFunction = app["instance-type"] == "carlos" ? readSBRPDataCarlos : readSBRPDataMatheus
  data, ids, data′, paths = readInstanceFunction(app)
  instance_name = split(basename(app["instance"]), ".")[1]
  function log(info)
    logColumns = ["model", "rootLP", "maxFlowCuts", "maxFlowCutsTime", "lazyCuts", "cost", "solverTime", "relativeGAP", "nodeCount"]
    println("instance", [" & :" * column for column in logColumns]...)
    println(instance_name, [" & " * (column in keys(info) ? string(info[column]) : "-") for column in logColumns]...)
  end
  println("|V| = $(length(data.D.V))")
  println("|A| = $(length(data.D.A))")
  println("|B| = $(length(data.B))")
  println("|V′| = $(length(data′.D.V))")
  println("|A′| = $(length(data′.D.A))")
  for b in data.B
    println("Block")
    for i in b
      println(data.D.V[ids[i]].pos_y, ", ", data.D.V[ids[i]].pos_x)
    end
  end
  # not solve
  app["nosolve"] && return
  # solve models
  #=
  # SBRP
  (model, x, info) = build_model_sbrp(data, app)
  optimize!(model)
  tour, info = gettour(data, x), merge(info, get_info(model))
  info["model"] = "SBRP"
  check_sbrp_sol(data, tour)
  log(info)
  if app["out"] != nothing
    writesol(app["out"], [ids[i] for i in tour if i != data.depot])
    writeGPX(app["out"] * ".gpx", [data.D.V[i] for i in tour if i != data.depot])
  end
  # SBRP complete 
  (model, x, y, info) = build_model_sbrp_complete(data′, app)
  optimize!(model)
  tour′, tour, info = gettour(data′, x), Array{Int64, 1}(), merge(info, get_info(model))
  check_sbrp_sol(data′, tour′)
  # replace compact paths
  for i in 2:length(tour′)
    push!(tour, tour′[i - 1])
    a = (tour′[i - 1], tour′[i])
    !in(a, data.D.A) && push!(tour, paths[a]...)
  end
  push!(tour, tour′[end])
  info["model"] = "SBRPComplete"
  check_sbrp_sol(data, tour)
#  println(tour)
  log(info)
  if app["out"] != nothing
    writesol(app["out"] * ".complete", [ids[i] for i in tour′ if i != data.depot])
    writeGPX(app["out"] * ".complete.gpx", [data.D.V[i] for i in tour′ if i != data.depot])
  end
  # TSP
  Vᵗ, Aᵗ, costs, Vb, Vb′ = build_atsp_instance(data′)
  (model, x, y) = build_model_atsp(Vᵗ, Aᵗ, costs)
  optimize!(model)
  println(objective_value(model))

  costs′ = from_atsp_to_tsp(Vᵗ, Aᵗ, costs)
  n = length(costs′)
  print("  ")
  [print(i, "    ") for i in floor(Int64, n/2):n]
  println()
  for i in 1:floor(Int64, n/2)
    print(i, " ")
    for j in floor(Int64, n/2):n
      d = costs′[i][j]
      log10 = Base.log(10, d)
      print(d, [" " for i in 1:(4 - (d != 0 ? floor(Int64, log10) : 0))]..., " ")
    end
    println()
  end

  print("   ")
  [print(i, "     ") for i in 1:floor(Int64, n/2)]
  println()
  for i in floor(Int64, n/2):n
    print(i, " ")
    for j in 1:floor(Int64, n/2)
      d = costs′[i][j]
      log10 = Base.log(10, d)
      print(d, [" " for i in 1:(4 - (d != 0 ? floor(Int64, log10) : 0))]..., " ")
    end
    println()
  end
  println(solve_tsp(costs′))
  
#    app["out"] != nothing && writesol(app["out"], data′, x, model, app)
  println("#######################SBRP_ATSP########################")
  V, A, costs, Vb, Vb′ = build_atsp_instance(data′)
  (model, x, y) = build_model_atsp(V, A, costs, data′.depot)
  if solve(model) 
    println(objective_value(model)) 
    tour = gettour(V, A, data′.depot, x)
    println(tour)
    for (i, b) in keys(Vb)
      println(i, " ", b, " ", Vb[(i, b)], " ", Vb′[(i, b)])
    end
    println(check_atsp_sol(tour, Vb, Vb′))
#    app["out"] != nothing && writesol(app["out"], data′, x, model, app)
  else
    println("Model infeasible or unknown")
  end
  =#
  (model, x, y, info) = build_model_sbrp_max(data′, app)
  optimize!(model)
  B = get_blocks(data, y)
  tour′, tour, info = gettour(data′, x, B), Array{Int64, 1}(), merge(info, get_info(model))
  [println(b) for b in B]
  check_sbrp_sol(data′, tour′, B)
  # replace compact paths
  for i in 2:length(tour′)
    push!(tour, tour′[i - 1])
    a = (tour′[i - 1], tour′[i])
    !in(a, data.D.A) && push!(tour, paths[a]...)
#    a in keys(paths) && push!(tour, paths[a]...)
  end
  push!(tour, tour′[end])
  info["model"] = "SBRPMax"
  check_sbrp_sol(data, tour, B)
  log(info)
  if app["out"] != nothing
    writesol(app["out"] * ".max", [ids[i] for i in tour if i != data.depot])
    writeGPX(app["out"] * ".max.gpx", [data.D.V[i] for i in tour if i != data.depot])
  end
  #=
  println("#######################SBRP MAX#########################")
  (model_max, x, y) = build_model_max_profit_sbrp(data, app)
  if solve(model_max)
    println(objective_value(model_max)) 
    app["out"] != nothing && writesol(app["out"] * "max", data, x, model)
  else
    println("Model infeasible or unknown")
  end
  println("########################################################")
  =#
end

function main(args)
   appfolder = dirname(@__FILE__)
   app = parse_commandline(args, appfolder)
   isnothing(app) && return
   if app["batch"] != nothing
      for line in readlines(app["batch"])
         if isempty(strip(line)) || strip(line)[1] == '#'
            continue
         end
         args_array = [String(s) for s in split(line)]
         app_line = parse_commandline(args_array, appfolder)
         run(app_line)
      end
   else
      run(app)
   end
end

#=
if isempty(ARGS)
   main(["--help"])
else
   main(ARGS)
end
=#

end
