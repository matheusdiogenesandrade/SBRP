using ArgParse

include("data.jl")
include("model.jl")
include("solve.jl")
include("sol.jl")

function parse_commandline(args_array::Array{String,1}, appfolder::String)
  s = ArgParseSettings(usage="  On interactive mode, call main([\"arg1\", ..., \"argn\"])", exit_after_help=false)
  @add_arg_table s begin
    "instance"
    help = "Instance file path"
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
  instance_name = split(basename(app["instance"]), ".")[1]
  data, ids, data′, paths = readSBRPData(app)
  println("|V| = $(length(data.D.V))")
  println("|A| = $(length(data.D.A))")
  println("|B| = $(length(data.B))")
  println("|V′| = $(length(data′.D.V))")
  println("|A′| = $(length(data′.D.A))")
  for b in data.B
    println(b)
  end
  # not solve
  app["nosolve"] && return
  # solve models
  println("#######################SBRP#############################")
  (model, x) = build_model_sbrp(data, app)
  if solve(model) 
    println(objective_value(model)) 
    tour = gettour(data, x)
    println(tour)
    check_sbrp_sol(data, tour)
#    app["out"] != nothing && writesol(app["out"], data, x, model, app)
  else
    println("Model infeasible or unknown")
  end
  println("#######################SBRP_COMPLETE####################")
  (model, x, y) = build_model_sbrp_complete(data′, app)
  if solve(model) 
    println(objective_value(model)) 
    tour = gettour(data′, x)
    println(tour)
    check_sbrp_sol(data′, tour)
#    app["out"] != nothing && writesol(app["out"], data′, x, model, app)
  else
    println("Model infeasible or unknown")
  end
  #=
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

if isempty(ARGS)
   main(["--help"])
else
   main(ARGS)
end
