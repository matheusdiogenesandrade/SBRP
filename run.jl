module Main

include("symbols.jl")
include("data.jl")
include("sol.jl")
include("nn_heuristic.jl")
include("model.jl")
include("solve.jl")
include("brkga.jl")
include("dd_sbrp.jl")

using .Data
using .Data.SBRP
using .Model
using .Model.ModelSBRP
using .Solution
using .NearestNeighborhoodHeuristic
using .BRKGA
using .DD_SBRP
using ArgParse

function parse_commandline(args_array::Array{String,1}, appfolder::String)
    s = ArgParseSettings(usage="  On interactive mode, call main([\"arg1\", ..., \"argn\"])", exit_after_help=false)
    @add_arg_table s begin
        "instance"
        help = "Instance file path"
        "--unitary-profits"
        help = "true if you want to consider the blocks' profitsas 1, and false otherwise"
        action = :store_true
        "--ip"
        help = "true if you want to run the I.P. model, and false otherwise"
        action = :store_true
        "--brkga"
        help = "true if you want to run the BRKGA, and false otherwise"
        action = :store_true
        "--brkga-conf"
        help = "BRKGA config file directory"
        default = "conf/config.conf"
        "--dd"
        help = "true if you want to run the Decision Diagram, and false otherwise"
        action = :store_true
        "--vehicle-time-limit"
        help = "Vehicle time limit in minutes"
        default = "120"
        "--instance-type"
        help = "Instance type (matheus|carlos)"
        default = "matheus"
        "--lb"
        help = "Lower bound"
        "--warm-start-solution"
        help = "Solution directory for Warm start"
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
        "--arcs-mtz"
        help = "MTZ on the arcs for the complete model"
        action = :store_true
        "--subcycle-separation"
        help = "Subcycle separation with max-flow at the B&B root node"
        action = :store_true
        "--y-integer"
        help = "Fix the variable y, for the complete model, when running the separation algorithm"
        action = :store_true
    end
    return parse_args(args_array, s)
end

# log function 
function log(info)
    logColumns = ["instance", "|V|", "|A|", "|B|", "T", "model", "initialLP", "yLP", "yLPTime", "intersectionCuts1", "intersectionCuts2", "intersectionCutsTime", "maxFlowLP", "maxFlowCuts", "maxFlowCutsTime", "lazyCuts", "cost", "solverTime", "relativeGAP", "nodeCount", "meters", "tourMinutes", "blocksMeters"]
    println([" & :" * column for column in logColumns]...)
    println([" & " * (column in keys(info) ? string(info[column]) : "-") for column in logColumns]...)
end

function ip(app::Dict{String, Any}, data::SBRPData, data′::SBRPData, paths::Dict{Tuple{Int, Int}, Vi})
    flush_println("###################SBRP####################")

    # create and solve model
    (model, x, y, info) = build_model_sbrp(data′, app)
    optimize!(model)

    # warm start checking
    if app["warm-start-solution"] != nothing && !has_values(model)

        # get solution
        tour′, B = Solution.readSolution(app["warm-start-solution"], data)

    else

	# get serviced blocks
	B = get_blocks(data, y) 

	# get solution for complete model
	tour′ = gettour(data′, x, B)

    end

    # check feasibility
    check_sbrp_sol(data′, tour′, B) 

    # get solution for original graph
#    tour = Vi()
#    for i in 2:length(tour′)
#        push!(tour, tour′[i - 1])
#        !in((tour′[i - 1], tour′[i]), data.D.A) && push!(tour, paths[(tour′[i - 1], tour′[i])]...)
#    end
#    push!(tour, tour′[end]) 

    # check feasibility
#    check_sbrp_sol(data, tour, B) 

    # log
    info = merge(info, get_info(model, data′, tour′, B), Dict{String, String}(
                                                                              "model" => "IP", 
                                                                              "instance" => app["instance_name"], 
                                                                              "|V|" => string(length(data′.D.V)), 
                                                                              "|A|" => string(length(data′.D.A)), 
                                                                              "|B|" => string(length(data.B)), 
                                                                              "T" => string(data.T))) 

    log(info) 

    # write solution
    if app["out"] != nothing
#        write_sol(app["out"], tour, data, B)
        write_sol(app["out"] * "_complete", tour′, data′, B)
        writeCostPerTime(COST_PER_TIME.list, app["out"] * "_cost_per_time")
    end

    flush_println("########################################################")
end

function brkga(app::Dict{String, Any}, data::SBRPData, data′::SBRPData, paths::Dict{Tuple{Int, Int}, Vi})
    flush_println("###################BRKGA####################")

    # solve model
    tour′, info, B, COST_PER_TIME = run_brkga(app["brkga-conf"], data′)

    # get tour 
    push!(tour′, data.depot)
    pushfirst!(tour′, data.depot)

    # check feasibility
    check_sbrp_sol(data′, tour′, B) 

    # get solution for original graph
#    tour = Vi()
#    for i in 2:length(tour′)
#        push!(tour, tour′[i - 1])
#        !in((tour′[i - 1], tour′[i]), data.D.A) && push!(tour, paths[(tour′[i - 1], tour′[i])]...)
#    end
#    push!(tour, tour′[end]) 

    # check feasibility
#    check_sbrp_sol(data, tour, B) 

    # log
    info = merge(info, Dict{String, String}(
                                            "model" => "BKRGA", 
                                            "instance" => app["instance_name"], 
                                            "|V|" => string(length(data′.D.V)), 
                                            "|A|" => string(length(data′.D.A)), 
                                            "|B|" => string(length(data.B)), 
                                            "T" => string(data.T),
#                                            "meters" => string(tour_distance(data, tour)),
#                                            "tourMinutes" => string(tour_time(data, tour, B)),
#                                            "blocksMeters" => string(sum(distance_block(data, block) for block in B))
                                            )) 
    log(info) 

    # write solution
    if app["out"] != nothing
#        write_sol(app["out"], tour, data, B)
        write_sol(app["out"] * "_complete", tour′, data′, B)
        writeCostPerTime(COST_PER_TIME, app["out"] * "_cost_per_time")
    end

    flush_println("########################################################")
end

function dd(app::Dict{String, Any}, data::SBRPData, data′::SBRPData, paths::Dict{Tuple{Int, Int}, Vi})
    flush_println("###################BRKGA####################")

    tour′, info, B = run_dd(data′)

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
                                            "model" => "DD", 
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
    if app["out"] != nothing
        write_sol(app["out"], tour, data, B)
        write_sol(app["out"] * "_complete", tour′, data′, B)
    end

    flush_println("########################################################")
end

function run(app::Dict{String,Any})
    flush_println("Application parameters:")
    [flush_println("  $arg  =>  $(repr(val))") for (arg,val) in app]

    # read instance
    instanceReader = app["instance-type"] == "matheus" ? readSBRPData : readSBRPDataCarlos

    data, data′, paths′ = instanceReader(app)
    app["instance_name"] = split(basename(app["instance"]), ".")[1]

    # instance data
    flush_println("|B| = $(length(data.B))")
    flush_println("|V| = $(length(data.D.V))")
    flush_println("|A| = $(length(data.D.A))")
    flush_println("|V′| = $(length(data′.D.V))")
    flush_println("|A′| = $(length(data′.D.A))")

    # set vehicle time limit
    data′.T = data.T = parse(Int, app["vehicle-time-limit"])

    # not solve
    app["nosolve"] && return

    # solve models
    app["ip"] && return ip(app, data, data′, paths′)
    app["brkga"] && return brkga(app, data, data′, paths′)
    app["dd"] && return dd(app, data, data′, paths′)

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
