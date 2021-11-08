using ArgParse

include("data.jl")
include("model.jl")
include("solve.jl")

function parse_commandline(args_array::Array{String,1}, appfolder::String)
   s = ArgParseSettings(usage="  On interactive mode, call main([\"arg1\", ..., \"argn\"])", exit_after_help=false)
   @add_arg_table s begin
      "instance"
         help = "Instance file path"
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
  data, ids = readSBRPData(app)
  # solve models
  println("#######################SBRP#############################")
  (model, x) = build_model_sbrp(data, app)
  if solve(model)
    println(objective_value(model)) 
    app["out"] != nothing && writesol(app["out"], data, x, model, app)
  else
    println("Model infeasible or unknown")
  end
  println("#######################SBRP MAX#########################")
  (model_max, x, y) = build_model_max_profit_sbrp(data, app)
  if solve(model_max)
    println(objective_value(model_max)) 
    app["out"] != nothing && writesol(app["out"] * "max", data, x, model)
  else
    println("Model infeasible or unknown")
  end
  println("########################################################")
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
