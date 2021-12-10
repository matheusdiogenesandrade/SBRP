module Main

include("data.jl")
include("model.jl")
include("solve.jl")
include("sol.jl")

using .Data
using .Data.SBRP
using .Model
using .Model.ModelSBRPMax
using .Model.ModelSBRPMaxComplete
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
  readInstanceFunction = app["instance-type"] == "carlos" ? readSBRPDataCarlos : readSBRPDataMatheus; data, ids, data′, paths = readInstanceFunction(app)
  instance_name = split(basename(app["instance"]), ".")[1]
  # log function 
  function log(info)
    logColumns = ["model", "rootLP", "maxFlowCuts", "maxFlowCutsTime", "lazyCuts", "cost", "solverTime", "relativeGAP", "nodeCount", "meters", "tourMinutes", "blocksMeters"]
    println("instance", [" & :" * column for column in logColumns]...)
    println(instance_name, [" & " * (column in keys(info) ? string(info[column]) : "-") for column in logColumns]...)
  end
  # instance data
  println("|V| = $(length(Set{Int64}(vcat([i for (i, j) in data.D.A], [j for (i, j) in data.D.A]))))")
  println("|A| = $(length(data.D.A))")
  println("|B| = $(length(data.B))")
  println("|V′| = $(length(Set{Int64}(vcat([i for (i, j) in data′.D.A], [j for (i, j) in data′.D.A]))))")
  println("|A′| = $(length(data′.D.A))")
  # not solve
  app["nosolve"] && return
  # solve models
  println("###################SBRP MAX#############################")
  (model, x, y, info) = build_model_sbrp_max(data, app); optimize!(model) # solve model
  B = get_blocks(data, y) # get serviced blocks
#  [println(block) for block in B]
  tour = gettour(data, x, B); check_sbrp_sol(data, tour, B) # get tour and check feasibility
  info = merge(info, get_info(model, data, tour, B), Dict{String, String}("model" => "SBRPMax")) # update info
  log(info) # log
  app["out"] != nothing && (writesol(app["out"] * ".max", [ids[i] for i in tour if i != data.depot]), writeGPX(app["out"] * ".max.gpx", [data.D.V[i] for i in tour if i != data.depot])) # write output files
  println("########################################################")
  #
  println("###################SBRP MAX Complete####################")
  (model, x, y, info) = build_model_sbrp_max_complete(data′, app); optimize!(model) # solve model
  B = get_blocks(data, y) # get serviced blocks
#  [println(block) for block in B]
  tour′ = gettour(data′, x, B); check_sbrp_sol(data′, tour′, B) # get tour and check feasibility
  info = merge(info, get_info(model, data′, tour′, B), Dict{String, String}("model" => "SBRPMaxComplete")) # update info
  tour = Array{Int64, 1}(); [(push!(tour, tour′[i - 1]); !in((tour′[i - 1], tour′[i]), data.D.A) && push!(tour, paths[(tour′[i - 1], tour′[i])]...)) for i in 2:length(tour′)]; push!(tour, tour′[end]) # replace compact paths
  check_sbrp_sol(data, tour, B) # check feasibility
  log(info) # log
  app["out"] != nothing && (writesol(app["out"] * ".complete.max", [ids[i] for i in tour if i != data.depot]), writeGPX(app["out"] * ".complete.max.gpx", [data.D.V[i] for i in tour if i != data.depot])) # write output files
  println("########################################################")
  #=
  i = 1
  while "iteration_" * string(i) * "_time" in keys(info)
    println(i, ": ", info["iteration_" * string(i) * "_time"], ", ", info["iteration_" * string(i) * "_cuts"])
    i += 1
  end
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
