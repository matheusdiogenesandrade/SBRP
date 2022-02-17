module Main

include("symbols.jl")
include("data.jl")
include("nn_heuristic.jl")
include("model.jl")
include("solve.jl")
include("sol.jl")
include("brkga.jl")

using .Data
using .Data.SBRP
using .Model
using .Model.ModelSBRPMax
using .Model.ModelSBRPMaxComplete
using .Solution
using .NearestNeighborhoodHeuristic
using .BRKGA
using ArgParse

function parse_commandline(args_array::Array{String,1}, appfolder::String)
  s = ArgParseSettings(usage="  On interactive mode, call main([\"arg1\", ..., \"argn\"])", exit_after_help=false)
  @add_arg_table s begin
    "instance"
    help = "Instance file path"
    "--instance-type"
    help = "Instance type (carlos, matheus)"
    default = "carlos"
    "--complete"
    help = "true if you want to run the model with the complete digraph, and false otherwise"
    action = :store_true
    "--no-one-degree-path"
    help = "true if you want to run the model with the no one degree path digraph, and false otherwise"
    action = :store_true
    "--normal"
    help = "true if you want to run the model with the original digraph, and false otherwise"
    action = :store_true
    "--brkga"
    help = "true if you want to run the BRKGA, and false otherwise"
    action = :store_true
    "--brkga-conf"
    help = "BRKGA config file directory"
    default = "conf/config.conf"
    "--vehicle-time-limit"
    help = "Vehicle time limit in minutes"
    default = "120"
    "--nosolve"
    help = "Not solve flag"
    action = :store_true
    "--out"
    help = "Path to write the solution found"
    "--batch"
    help = "Batch file path"
    "--intersection-cuts"
    help = "Intersection model for the complete model"
    action = :store_true
    "--y-integer"
    help = "Fix the variable y, for the complete model, when running the separation algorithm"
    action = :store_true
  end
  return parse_args(args_array, s)
end

# log function 
function log(info)
  logColumns = ["instance", "|V|", "|A|", "|B|", "T", "model", "initialLP", "yLP", "yLPTime", "intersectionCuts", "intersectionCutsTime", "maxFlowLP", "maxFlowCuts", "maxFlowCutsTime", "lazyCuts", "cost", "solverTime", "relativeGAP", "nodeCount", "meters", "tourMinutes", "blocksMeters"]
  println([" & :" * column for column in logColumns]...)
  println([" & " * (column in keys(info) ? string(info[column]) : "-") for column in logColumns]...)
end

function write_sol(app, tour, data)
  if app["out"] != nothing 
    writesol(app["out"] * ".sol", [i for i in tour if i != data.depot])
    # write output files
    writeGPX(app["out"] * ".gpx", [data.D.V[i] for i in tour if i != data.depot]) 
  end
end

function sbrp_max(app::Dict{String, Any}, data::SBRPData)
  println("###################SBRP MAX#############################")

  # create model
  (model, x, y, z, info) = build_model_sbrp_max(data, app)

  # solve model
  optimize!(model)
  println(objective_value(model))

  # get serviced blocks
  B = get_blocks(data, y)
#  [println(block) for block in B]

  # get tour
  tour = gettour(data, x, B)

  # check feasibility
  check_sbrp_sol(data, tour, B)

  # update info
  info = merge(info, get_info(model, data, tour, B), Dict{String, String}(
              "model" => "SBRPMax", 
              "instance" => app["instance_name"], 
              "|V|" => string(length(Set{Int}(vcat([i for (i, j) in data.D.A], [j for (i, j) in data.D.A])))), 
              "|A|" => string(length(data.D.A)), 
              "|B|" => string(length(data.B)), 
              "T" => string(data.T))) 

  # log
  log(info)

  # solution
  write_sol(app, tour, data)

  println("########################################################")
  #=
  i = 1
  while "iteration_" * string(i) * "_time" in keys(info)
    println(i, ": ", info["iteration_" * string(i) * "_time"], ", ", info["iteration_" * string(i) * "_cuts"])
    i += 1
  end
  =#
end
 
function sbrp_max_complete(app::Dict{String, Any}, data::SBRPData, data′::SBRPData, paths::Dict{Tuple{Int, Int}, Vi})
  flush_println("###################SBRP MAX Complete####################")

  # create and solve model
  (model, x, y, info) = build_model_sbrp_max_complete(data′, app)
  optimize!(model) 

  # get serviced blocks
  B = get_blocks(data, y) 
  
  # get solution for complete model
  tour′ = gettour(data′, x, B)

  # check feasibility
  check_sbrp_sol(data′, tour′, B) 

  # get solution for original graph
  tour = Vi()
  for i in 2:length(tour′)
    push!(tour, tour′[i - 1])
    !in((tour′[i - 1], tour′[i]), data.D.A) && push!(tour, paths[(tour′[i - 1], tour′[i])]...)
  end
  push!(tour, tour′[end]) 

  # check feasibility
  check_sbrp_sol(data, tour, B) 

  # log
  info = merge(info, get_info(model, data′, tour′, B), Dict{String, String}(
               "model" => "SBRPMaxComplete", 
               "instance" => app["instance_name"], 
               "|V|" => string(length(data′.D.V)), 
               "|A|" => string(length(data′.D.A)), 
               "|B|" => string(length(data.B)), 
               "T" => string(data.T))) 

  log(info) 

  # write solution
  write_sol(app, tour, data)

  flush_println("########################################################")
  #=
  i = 1
  while "iteration_" * string(i) * "_time" in keys(info)
    println(i, ": ", info["iteration_" * string(i) * "_time"], ", ", info["iteration_" * string(i) * "_cuts"])
    i += 1
  end
  =#
end

function brkga(app::Dict{String, Any}, data::SBRPData, data′::SBRPData, paths::Dict{Tuple{Int, Int}, Vi})
  flush_println("###################BRKGA####################")

  # solve model
  tour′, info, B = run_brkga(app["brkga-conf"], data′)

  # get tour 
  push!(tour′, data.depot)
  pushfirst!(tour′, data.depot)

  # check feasibility
  check_sbrp_sol(data′, tour′, B) 

  # get solution for original graph
  tour = Vi()
  for i in 2:length(tour′)
    push!(tour, tour′[i - 1])
    !in((tour′[i - 1], tour′[i]), data.D.A) && push!(tour, paths[(tour′[i - 1], tour′[i])]...)
  end
  push!(tour, tour′[end]) 

  # check feasibility
  check_sbrp_sol(data, tour, B) 

  # log
  info = merge(info, Dict{String, String}(
                                          "model" => "BKRGA", 
                                          "instance" => app["instance_name"], 
                                          "|V|" => string(length(data′.D.V)), 
                                          "|A|" => string(length(data′.D.A)), 
                                          "|B|" => string(length(data.B)), 
                                          "T" => string(data.T),
                                          "meters" => string(tour_distance(data, tour)),
                                          "tourMinutes" => string(tour_time(data, tour, B)),
                                          "blocksMeters" => string(sum(distance_block(data, block) for block in B))
                                         )) 
  log(info) 

  # write solution
  write_sol(app, tour, data)

  flush_println("########################################################")
end

function sbrp_max_no_one_degree_path(app::Dict{String, Any}, data::SBRPData, data″::SBRPData, paths″::Dict{Arc, Vi})

  flush_println("#############SBRP MAX NO ONE DEGREE PATH################")

  # create model
  (model, x, y, z, info) = build_model_sbrp_max(data″, app)

  # solve model
  optimize!(model) 
  println(objective_value(model))

  # get serviced blocks
  B = get_blocks(data, y) 
#  [println(block) for block in B]

  # get tour for the processed digraph
  tour″ = gettour(data″, x, B)

  # check feasibility
  check_sbrp_sol(data″, tour″, B) 

  # get solution for the original graph
  tour = Vi()
  for i in 2:length(tour″)
    a = curr, next = tour″[i - 1], tour″[i]
    # store node
    push!(tour, curr)

    # check if there is some path
    !in(a, data.D.A) && push!(tour, paths″[a]...)
  end
  push!(tour, tour″[end]) 

  # check feasibility
  check_sbrp_sol(data, tour, B) 

  # log
  info = merge(info, get_info(model, data″, tour″, B), Dict{String, String}(
               "model" => "SBRPMaxNoOneDegreePath", 
               "instance" => app["instance_name"], 
               "|V|" => string(length(data″.D.V)), 
               "|A|" => string(length(data″.D.A)), 
               "|B|" => string(length(data″.B)), 
               "T" => string(data″.T))) 

  log(info) 

  # write solution
  write_sol(app, tour, data)

  flush_println("########################################################")
  #=
  i = 1
  while "iteration_" * string(i) * "_time" in keys(info)
    println(i, ": ", info["iteration_" * string(i) * "_time"], ", ", info["iteration_" * string(i) * "_cuts"])
    i += 1
  end
  =#
end

function run(app::Dict{String,Any})
  flush_println("Application parameters:")
  [flush_println("  $arg  =>  $(repr(val))") for (arg,val) in app]

  # read instance
  readInstanceFunction = app["instance-type"] == "carlos" ? readSBRPDataCarlos : readSBRPDataMatheus
  data, data′, paths′, data″, paths″ = readInstanceFunction(app)
  app["instance_name"] = split(basename(app["instance"]), ".")[1]

  # instance data
  flush_println("|B| = $(length(data.B))")
  flush_println("|V| = $(length(data.D.V))")
  flush_println("|A| = $(length(data.D.A))")
  flush_println("|V′| = $(length(data′.D.V))")
  flush_println("|A′| = $(length(data′.D.A))")
  flush_println("|V″| = $(length(data″.D.V))")
  flush_println("|A″| = $(length(data″.D.A))")

  # set vehicle time limit
  data″.T = data′.T = data.T = parse(Int, app["vehicle-time-limit"])

  # not solve
  app["nosolve"] && return

  # solve models
  app["normal"] && sbrp_max(app, data)
  app["no-one-degree-path"] && sbrp_max_no_one_degree_path(app, data, data″, paths″)
  app["complete"] && sbrp_max_complete(app, data, data′, paths′)
  app["brkga"] && brkga(app, data, data′, paths′)

end

end

using .Main

function main(args)
  appfolder = dirname(@__FILE__)
  app = Main.parse_commandline(args, appfolder)
  isnothing(app) && return
  if app["batch"] != nothing 
    for line in readlines(app["batch"]) 
      if !isempty(strip(line)) && strip(line)[1] != '#'
        Main.run(Main.parse_commandline([String(s) for s in split(line)], appfolder))
      end
    end
  else
    Main.run(app)
  end
end

if isempty(ARGS)
   main(["--help"])
else
   main(ARGS)
end
