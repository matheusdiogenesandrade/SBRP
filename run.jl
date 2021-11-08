using ArgParse

include("data.jl")
include("model.jl")

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

   instance_name = split(basename(app["instance"]), ".")[1]

   data = readSBRPData(app)

   # (model, x, y) = build_model(data, app)
   # solve model
   # (status, solution_found) = optimize!(optimizer)
   solution_found = false
   if solution_found
      println("########################################################")
      #=
      (app["out"] != nothing) && write(f, "Cost: $(get_objective_value(optimizer))\n")
      (app["out"] != nothing) && close(f)
      =#
      println("########################################################")
   elseif status == :Infeasible
      println("Problem infeasible")
   else
      println("Solution not found")
   end
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
