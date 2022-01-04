module Main


include("symbols.jl")
include("data.jl")
include("nn_heuristic.jl")
include("model.jl")
include("solve.jl")
include("sol.jl")

using .Symbols
using .Data
using .Data.SBRP
using .Model
using .Model.ModelSBRPMax
using .Model.ModelSBRPMaxComplete
using .Solution
using .NearestNeighborhoodHeuristic
using ArgParse

function parse_commandline(args_array::Array{String,1}, appfolder::String)
  s = ArgParseSettings(usage="  On interactive mode, call main([\"arg1\", ..., \"argn\"])", exit_after_help=false)
  @add_arg_table s begin
    "instance"
    help = "Instance file path"
    "--instance-type"
    help = "Instance type (carlos, matheus)"
    default = "carlos"
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
function sbrp_max(app::Dict{String, Any}, data::SBRPData, ids::Dict{Int64, Int64})
  println("###################SBRP MAX#############################")
  (model, x, y, info) = build_model_sbrp_max(data, app); optimize!(model) # solve model
  B = get_blocks(data, y) # get serviced blocks
#  [println(block) for block in B]
  tour = gettour(data, x, B); check_sbrp_sol(data, tour, B) # get tour and check feasibility
  info = merge(info, get_info(model, data, tour, B), Dict{String, String}("model" => "SBRPMax", "instance" => app["instance_name"], "|V|" => string(length(Set{Int64}(vcat([i for (i, j) in data.D.A], [j for (i, j) in data.D.A])))), "|A|" => string(length(data.D.A)), "|B|" => string(length(data.B)), "T" => string(data.T))) # update info
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
 
function sbrp_max_complete(app::Dict{String, Any}, data::SBRPData, data′::SBRPData, ids::Dict{Int64, Int64}, paths::Dict{Tuple{Int64, Int64}, Array{Int64}})
  println("###################SBRP MAX Complete####################")
  flush(stdout)
  (model, x, y, info) = build_model_sbrp_max_complete(data′, app); optimize!(model) # solve model
  B = get_blocks(data, y) # get serviced blocks
  #  [println(block) for block in B]
  tour′ = gettour(data′, x, B); check_sbrp_sol(data′, tour′, B) # get tour and check feasibility
  info = merge(info, get_info(model, data′, tour′, B), Dict{String, String}("model" => "SBRPMaxComplete", "instance" => app["instance_name"], "|V|" => string(length(Set{Int64}(vcat([i for (i, j) in data′.D.A], [j for (i, j) in data′.D.A])))), "|A|" => string(length(data′.D.A)), "|B|" => string(length(data.B)), "T" => string(data.T))) # update info
  tour = Array{Int64, 1}(); [(push!(tour, tour′[i - 1]); !in((tour′[i - 1], tour′[i]), data.D.A) && push!(tour, paths[(tour′[i - 1], tour′[i])]...)) for i in 2:length(tour′)]; push!(tour, tour′[end]) # replace compact paths
  check_sbrp_sol(data, tour, B) # check feasibility
  log(info) # log
  write_sol(app, tour, ids, data)
  println("########################################################")
  flush(stdout)
  #=
  i = 1
  while "iteration_" * string(i) * "_time" in keys(info)
    println(i, ": ", info["iteration_" * string(i) * "_time"], ", ", info["iteration_" * string(i) * "_cuts"])
    i += 1
  end
  =#
end

function run(app::Dict{String,Any})
  println("Application parameters:")
  [println("  $arg  =>  $(repr(val))") for (arg,val) in app]
  flush(stdout)
  # read instance
  readInstanceFunction = app["instance-type"] == "carlos" ? readSBRPDataCarlos : readSBRPDataMatheus; data, ids, data′, paths = readInstanceFunction(app)
  app["instance_name"] = split(basename(app["instance"]), ".")[1]
  # instance data
  println("|V| = $(length(Set{Int64}(vcat([i for (i, j) in data.D.A], [j for (i, j) in data.D.A]))))")
  println("|A| = $(length(data.D.A))")
  println("|B| = $(length(data.B))")
  println("|V′| = $(length(Set{Int64}(vcat([i for (i, j) in data′.D.A], [j for (i, j) in data′.D.A]))))")
  println("|A′| = $(length(data′.D.A))")
  flush(stdout)
  # set vehicle time limit
  data′.T = data.T = parse(Int64, app["vehicle-time-limit"])
  # not solve
  app["nosolve"] && return
  # solve models
#  sbrp_max(app, data, ids)
  sbrp_max_complete(app, data, data′, ids, paths)
  flush(stdout)
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
