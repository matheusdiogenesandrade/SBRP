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
    help = "true if you want to run the model with complete digraph, and false otherwise"
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

function write_sol(app, tour, ids, data)
  app["out"] != nothing && (writesol(app["out"] * ".max", [ids[i] for i in tour if i != data.depot]), writeGPX(app["out"] * ".max.gpx", [data.D.V[i] for i in tour if i != data.depot])) # write output files
end

#=
function sbrp_max(app::Dict{String, Any}, data::SBRPData, ids::Dict{Int, Int})
  println("###################SBRP MAX#############################")
  (model, x, y, info) = build_model_sbrp_max(data, app); optimize!(model) # solve model
  println(objective_value(model))
  B = get_blocks(data, y) # get serviced blocks
#  [println(block) for block in B]
  tour = gettour(data, x, B); check_sbrp_sol(data, tour, B) # get tour and check feasibility
  info = merge(info, get_info(model, data, tour, B), Dict{String, String}("model" => "SBRPMax", "instance" => app["instance_name"], "|V|" => string(length(Set{Int}(vcat([i for (i, j) in data.D.A], [j for (i, j) in data.D.A])))), "|A|" => string(length(data.D.A)), "|B|" => string(length(data.B)), "T" => string(data.T))) # update info
  log(info) # log
  write_sol(app, tour, ids, data)
  println("########################################################")
  #=
  i = 1
  while "iteration_" * string(i) * "_time" in keys(info)
    println(i, ": ", info["iteration_" * string(i) * "_time"], ", ", info["iteration_" * string(i) * "_cuts"])
    i += 1
  end
  =#
end
=#
 
function sbrp_max_complete(app::Dict{String, Any}, data::SBRPData, data′::SBRPData, ids::Dict{Int, Int}, paths::Dict{Tuple{Int, Int}, Vi})
  flush_println("###################SBRP MAX Complete####################")

  # create and solve model
  (model, x, y, info) = build_model_sbrp_max_complete(data′, app); 
  optimize!(model) 

  # get serviced blocks
  B = get_blocks(data, y) 
  #[flush_println(block) for block in B]

  # get solution for complete model
  tour′ = gettour(data′, x, B); check_sbrp_sol(data′, tour′, B) # get tour and check feasibility

  # get solution for original graph
  tour = Vi()
  [(push!(tour, tour′[i - 1]); !in((tour′[i - 1], tour′[i]), data.D.A) && push!(tour, paths[(tour′[i - 1], tour′[i])]...)) for i in 2:length(tour′)]
  push!(tour, tour′[end]) 

  # check feasibility
  check_sbrp_sol(data, tour, B) 

  # log
  info = merge(info, get_info(model, data′, tour′, B), Dict{String, String}(
               "model" => "SBRPMaxComplete", 
               "instance" => app["instance_name"], 
               "|V|" => string(length(Set{Int}(vcat([i for (i, j) in data′.D.A], [j for (i, j) in data′.D.A])))), 
               "|A|" => string(length(data′.D.A)), 
               "|B|" => string(length(data.B)), 
               "T" => string(data.T))) 
  log(info) 

  # write solution
  write_sol(app, tour, ids, data)

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
  data, ids, data′, paths = readInstanceFunction(app)
  app["instance_name"] = split(basename(app["instance"]), ".")[1]
  # instance data
  flush_println("|V| = $(length(Set{Int}(vcat([i for (i, j) in data.D.A], [j for (i, j) in data.D.A]))))")
  flush_println("|A| = $(length(data.D.A))")
  flush_println("|B| = $(length(data.B))")
  flush_println("|V′| = $(length(Set{Int}(vcat([i for (i, j) in data′.D.A], [j for (i, j) in data′.D.A]))))")
  flush_println("|A′| = $(length(data′.D.A))")
  flush(stdout)
  # set vehicle time limit
  data′.T = data.T = parse(Int, app["vehicle-time-limit"])
  # not solve
  app["nosolve"] && return
  # solve models
#  sbrp_max(app, data, ids)
  app["complete"] && sbrp_max_complete(app, data, data′, ids, paths)
  app["brkga"] && run_brkga(app, data′)
end

end

using .Main

function main(args)
  appfolder = dirname(@__FILE__)
  app = Main.parse_commandline(args, appfolder)
  isnothing(app) && return
  app["batch"] != nothing &&  ([Main.run(Main.parse_commandline([String(s) for s in split(line)], appfolder)) for line in readlines(app["batch"]) if !isempty(strip(line)) && strip(line)[1] != '#']; return)
  Main.run(app)
end

if isempty(ARGS)
   main(["--help"])
else
   main(ARGS)
end
