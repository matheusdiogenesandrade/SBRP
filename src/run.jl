include("symbols.jl")
include("data.jl")
include("sol.jl")
include("nn_heuristic.jl")
include("model.jl")
include("brkga.jl")

using ArgParse
using Logging

# log
const LOG_FILE = open("logs/log", "w+")
global_logger(ConsoleLogger(LOG_FILE, Logging.Debug))
#disable_logging(Debug) 

function parse_commandline(args_array::Vector{String}, appfolder::String)::Union{Nothing, Dict{String, Any}}
    s::ArgParseSettings = ArgParseSettings(usage="  On interactive mode, call main([\"arg1\", ..., \"argn\"])", exit_after_help=false)

    @add_arg_table s begin
        "instance"
        help = "Instance file path"
        "--cluster-size-profits"
        help = "true if you want to consider the blocks' profits as the size of the block, and false otherwise"
        action = :store_true
        "--unitary-profits"
        help = "true if you want to consider the blocks' profits as 1, and false otherwise"
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
        "--vehicle-time-limit"
        help = "Vehicle time limit in minutes"
        default = "120"
        "--instance-type"
        help = "Instance type (matheus|carlos)"
        default = "matheus"
        "--nosolve"
        help = "Not solve flag"
        action = :store_true
        "--out"
        help = "Path to write the solution found"
        "--batch"
        help = "Batch file path"
        "--intersection-cuts"
        help = "Intersection cuts for the complete model"
        action = :store_true
        "--subcycle-separation"
        help = "Strategy (first|best|all|none)"
        help = "Subcycle separation with max-flow at the B&B root node"
        default = "all"
        "--y-integer"
        help = "Fix the variable y, for the complete model, when running the separation algorithm"
        action = :store_true
        "--z-integer"
        help = "Fix the variable z, for the complete model, when running the separation algorithm"
        action = :store_true
        "--w-integer"
        help = "Fix the variable w, for the complete model, when running the separation algorithm"
        action = :store_true
    end

    return parse_args(args_array, s)
end

# log function 
function log(app::Dict{String, Any}, info::Dict{String, String})

    columns::Vector{String} = ["instance", "|V|", "|A|", "|B|", "T", "model", "initialLP", "yLP", "yLPTime", "zLP", "zLPTime", "wLP", "wLPTime" , "maxFlowLP", "maxFlowCuts", "maxFlowCutsTime", "lazyCuts", "cost", "solverTime", "relativeGAP", "nodeCount", "integerCount", "meters", "tourMinutes", "blocksMeters", "numVisitedBlocks", "intersectionCutsTime", "intersectionCuts1", "intersectionCuts2"]

    info["instance"] = last(split(app["instance"], "/"; keepempty = false))
    info["instance"] = first(split(info["instance"], "."; keepempty = false))

    selected_columns = filter(column -> column in keys(info), columns)

    # log
    @info join(selected_columns, ",")
    @info join(map(column -> info[column], selected_columns), ",")

    # console dump
    println(join(selected_columns, ","))
    println(join(map(column -> info[column], selected_columns), ","))
    flush(stdout)

end

#=
Get solution for original graph
    input: 
        - tour::Vector{Integer} is the solution route
        - A::Vector{Pair{Int, Int}} is the list of arcs `tour` was built on
        - paths::Dict{Pair{Int, Int}, Vector{Int}} is the relation of paths, from the original digraph, between two nodes belonging to the blocks

    output:
        - ori_tour::Vector{Int} is the solution route using arcs from the original digraph
=# 
function retrieveOriginalDigraphSolution(
        tour::Vi,
        A::Arcs,
        paths::Dict{Arc, Vi}
    )::Vi

    # optimization
    A_set::ArcsSet = ArcsSet(A)

    # original route
    ori_tour::Vi = vcat(
         map(
             (i, j)::Arc -> (i, j) in A_set ? vcat(i, paths[(i, j)]) : [], 
             zip(tour[begin:end - 1], tour[2:end])
            )
        )
    push!(ori_tour, tour[end]) 

    # return
    return ori_tour
end

function getPathsForIntersectionCuts(data::SBRPData, info::Dict{String, String})::Tuple{VVi, ArcsSet}

    # get
    elapsed_time::Float64 = @elapsed maximal_paths::VVi, invalid_arcs::ArcsSet = getMaximalPathsAndInvalidArcs(data) 

    # paths statistics
    paths_lengths::Vi = map(maximal_path::Vi -> length(maximal_path), maximal_paths)

    num_paths::Int   = length(maximal_paths)
    average::Float64 = reduce(+, paths_lengths) / num_paths
    max_length::Int  = maximum(paths_lengths)
    min_length::Int  = minimum(paths_lengths)
    std::Float64     = reduce(+, map(path_length::Int -> abs(path_length - average), paths_lengths)) / num_paths

    # store
    info["#PreprocArcs"]    = string(length(invalid_arcs))
    info["#IntersectPaths"] = string(num_paths)
    info["AVGIntersectLen"] = string(average)
    info["MaxIntersectLen"] = string(max_length)
    info["MinIntersectLen"] = string(min_length)
    info["STDIntersectLen"] = string(std)
    info["IntersecTime"]    = string(elapsed_time)

    # return
    return maximal_paths, invalid_arcs
end

#=
Run complete digraph IP formulation for the SBRP
    input: 
        - app::Dict{String, Any} is the command line arguments relation
        - data::SBRPData is the SBRP instance built on a complete digraph
=# 
function completeDigraphIPModel(app::Dict{String, Any}, data::SBRPData)

    @info "###################SBRP####################"

    # intersection cuts
    info::Dict{String, String} = Dict{String, String}()
    maximal_paths::VVi         = VVi()

    # create and solve model
    solution::SBRPSolution, info_::Dict{String, String} = runCOPCompleteDigraphIPModel(data, app)

    # merge
    merge!(info, info_)

    # check feasibility
    checkSBRPSolution(data, solution) 

    # log
    info["model"]        = "IP"
    info["|V|"]          = string(length(data.D.V))
    info["|A|"]          = string(length(data.D.A))
    info["|B|"]          = string(length(data.B))
    info["T"]            = string(data.T)
    info["meters"]       = string(tourDistance(data, solution.tour))
    info["tourMinutes"]  = string(tourTime(data, solution))
    info["blocksMeters"] = string(sum(map(block::Vi -> blockDistance(data, block), solution.B)))
    info["numVisitedBlocks"] = string(length(solution.B))

    log(app, info) 

    # write solution
    solution_dir::Union{String, Nothing} = app["out"]

    if solution_dir != nothing
        writeSolution(solution_dir * "_cop_ip_model", data, solution)
    end

    @info "########################################################"
end

#=
Run BRKGA algorithm
    input: 
        - app::Dict{String, Any} is the command line arguments relation
        - data::SBRPData is the SBRP instance built on a complete digraph
=# 
function BRKGAModel(app::Dict{String, Any}, data::SBRPData)

    @info "###################BRKGA####################"

    # solve model
    solution::SBRPSolution, info::Dict{String, String} = runCOPBRKGAModel(data, app)

    # check feasibility
    checkSBRPSolution(data, solution) 

    # log
    info["model"]        = "BRKGA"
    info["|V|"]          = string(length(data.D.V))
    info["|A|"]          = string(length(data.D.A))
    info["|B|"]          = string(length(data.B))
    info["T"]            = string(data.T)
    info["meters"]       = string(tourDistance(data, solution.tour))
    info["tourMinutes"]  = string(tourTime(data, solution))
    info["blocksMeters"] = string(sum(map(block::Vi -> blockDistance(data, block), solution.B)))
    info["numVisitedBlocks"] = string(length(solution.B))

    log(app, info) 

    # write solution
    solution_dir::Union{String, Nothing} = app["out"]
    if solution_dir != nothing
        writeSolution(solution_dir * "_brkga", data, solution)
    end

    @info "########################################################"
end

function run(app::Dict{String,Any})
    @info "Application parameters:"

    for (arg, val) in app
        @info "  $arg  =>  $(repr(val))"
    end

    # read instance
    data::Union{SBRPData, Nothing} = nothing

    if app["instance-type"] == "matheus"
        data = readSBRPData(app)
    elseif app["instance-type"] == "carlos"
        data = readSBRPDataCarlos(app)
    else
        error(@sprintf("Invalid instance type: %s", app["instance-type"]))
    end

    # instance data
    @info "|B| = $(length(data.B))"
    @info "|V| = $(length(data.D.V))"
    @info "|A| = $(length(data.D.A))"

    # set vehicle time limit
    data.T = parse(Int, app["vehicle-time-limit"])

    # not solve
    if app["nosolve"] 

        info::Dict{String, String} = Dict{String, String}()

        # log
        info["|V|"] = string(length(data.D.V))
        info["|A|"] = string(length(data.D.A))
        info["|B|"] = string(length(data.B))
        info["T"]   = string(data.T)

        log(app, info)
    elseif app["ip"]
        completeDigraphIPModel(app, data)
    elseif app["brkga"] 
        BRKGAModel(app, data)
    end
end

function main(args)
    appfolder::String = dirname(@__FILE__)

    app::Union{Nothing, Dict{String, Any}} = parse_commandline(args, appfolder)

    isnothing(app) && return

    if app["batch"] != nothing 
        for line in readlines(app["batch"]) 
            if !isempty(strip(line)) && strip(line)[1] != '#'

                run(parse_commandline(map(s::Any -> String(s), split(line)), appfolder))
                GC.gc()  # Force the GC to run

            end
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
