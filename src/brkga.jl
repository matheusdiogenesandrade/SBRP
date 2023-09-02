using BrkgaMpIpr
using ConfParser
using Dates
using Printf

"""
    @enum StopRule

Controls stop criteria. Stops either when:
- a given number of `GENERATIONS` is given;
- or a `TARGET` value is found;
- or no `IMPROVEMENT` is found in a given number of iterations.
"""
@enum StopRule begin
    GENERATIONS = 0
    TARGET = 1
    IMPROVEMENT = 2
end

"""
    function parse(::Type{StopRule}, value::String)::StopRule

Parse `value` into a `StopRule`.
"""
function parseRule(::Type{StopRule}, value::String)::StopRule

    @debug "Parsing GA rules"

    local_value::Char = uppercase(strip(value)[1])

    dict::Dict{Char, StopRule} = Dict{Char, StopRule}('G' => GENERATIONS, 'T' => TARGET, 'I' => IMPROVEMENT)

    local_value in keys(dict) && return dict[local_value]

    throw(ArgumentError("cannot parse $value as StopRule"))
end

function dijkstra(data::SBRPData, idxs_blocks::Vi)::Tuple{Float64, Dict{Tuple{Int, Int}, Float64}}

    @debug "Dijkstra algorithm"


    # attrs
    
    @debug "Getting params"

    m::Int = length(idxs_blocks)
    B::VVi         = data.B
    A::Arcs        = data.D.A
    T::Float64     = data.T
    first_idx::Int = first(idxs_blocks) 

    # dijkstra distance <<idx_block, node>, min path time> stores the min path time from the depot to the node `node' in the block at `idx_block'
    consumed_times::Dict{Tuple{Int, Int}, Float64} = Dict{Tuple{Int, Int}, Float64}(reduce(vcat,
                                                                                           map(idx_block::Int -> 
                                                                                               map(i::Int -> (idx_block, i)::Tuple{Int, Int} => T + 1, B[idx_block]),
                                                                                               idxs_blocks
                                                                                              )
                                                                                          ))

    # initial times

    @debug "Initial times"

    first_block_time::Float64 = blockTime(data, B[first_idx])
    for i::Int in B[first_idx]
        consumed_times[(first_idx, i)] = first_block_time
    end

    # BFS

    @debug "BFS"

    for i::Int in 1:m - 1

        # get indexes

        @debug "Get indexes"

        idx_curr_block::Int = idxs_blocks[i]
        idx_next_block::Int = idxs_blocks[i + 1] 

        # get curr and next blocks

        @debug "Get curr and next block"

        curr_block::Vi = B[idx_curr_block]
        next_block::Vi = B[idx_next_block]

        # calculate next block service time
        
        @debug "Get next block service time"

        next_block_time::Float64 = blockTime(data, next_block)

        # intersecting nodes
       
        @debug "Get intersecting nodes"

        intersection::Si = Si(intersect(curr_block, next_block))

        for j::Int in intersection
            consumed_times[(idx_next_block, j)] = next_block_time + consumed_times[(idx_curr_block, j)]
        end

        # non intersecting nodes
        
        @debug "Get non-intersecting nodes"

        for i::Int in curr_block
            for j::Int in next_block
                if Arc(i, j) in keys(data.D.distance)

                    consumed_times[(idx_next_block, j)] = min(consumed_times[(idx_next_block, j)], consumed_times[(idx_curr_block, i)] + time(data, Arc(i, j)) + next_block_time)

                end
            end
        end
    end

    @debug "Get last block index"

    last_block_idx::Int = idxs_blocks[findlast(idx::Int -> any(i::Int -> consumed_times[(idx, i)] <= T, B[idx]), idxs_blocks)]

    @debug "Return"

    candidates::Vi = filter(i::Int -> consumed_times[(last_block_idx, i)] <= T, B[last_block_idx])

    return minimum(map(i::Int -> consumed_times[(last_block_idx, i)], candidates)), consumed_times
end

function getLongestProfitRoute(data::SBRPData, consumed_times::Dict{Tuple{Int, Int}, Real}, idxs_blocks::Vi)

    #=
    # params
    V = data.D.V
    depot = data.depot
    B = data.B

    # blocks order
    blocks_order::VVi = map(idx::Int -> B[idx], idxs_blocks)
    pushfirst!(blocks_order, [depot])

    # build DAG
    # build nodes
    dag_V::Dict{Int, Vertex} = Dict{Int, Vertex}()

    block_node_to_id::Dict{Tuple{Vi, Int}, Int} = Dict{Tuple{Vi, Int}, Int}()
    id_to_block_node::Dict{Int, Tuple{Vi, Int}} = Dict{Int, Tuple{Vi, Int}}()

    id::Int = 0
    for block in blocks_order
    for i in block
    dat_V[id] = V[i]
    block_node_to_id[(block, i)] = id
    id_to_block_node[id] = (block, i)

    id += 1
    end
    end

    # build arcs
    dag_A::Arcs = Arcs()
    enumeration = enumerate(blocks_order)
    for idx_prev, block_prev in enumeration[1:end - 1]
    for _, block_next in enumeration[idx_prev + 1:end]
    for i in block_prev
    for j in block_next
    node_prev = block_node_to_id[(block_prev, i)]
    node_next = block_node_to_id[(block_next, j)]

    push!(dag_arcs, Arc(node_prev, node_next))
    end
    end
    end
    end

    # build revAdjList
    revAdjList::Dict{Int, Vi} = Dict{Int, Vi}(id => Vi())
    for (id_i, id_j) in dag_A 
    push!(revAdjList[id_j], id_i)
    end

    # times and profits
    times::ArcCostMap = ArcCostMap()
    profits::ArcCostMap = ArcCostMap()
    for (id_i, id_j) in dag_A
    i, block_i = id_to_block_node[id_i]
    j, block_j = id_to_block_node[id_j]

    times[(id_i, id_j)] = time(data, (i, j)) + time_block(data, block_j)
    profits[(id_i, id_j)] = data.profits[block_j]
    end

    # get topological sorting
    topological_order::Vi = Vi()
    for block in blocks_order
    for i in block
    id_i = block_node_to_id[(block, i)]

    push!(topological_order, id_i)
    end
    end

    #
    path::Vi = Vi()
    pred::Dict{Int, Int} = Dict{Int, Int}(i => i for i in keys(dag_V))
    incurr_time::Vector{Float64} = Vector{Float64}(zeros(length(pred)))
    incurr_profit::Vector{Float64} = Vector{Float64}(zeros(length(pred)))

    # get distances
    for v in topological_order
    for u in revAdjList[v]

    newdist = dist[u] + times[(u, v)]
    newprofit = 

    if newdist > dist[v]
    dist[v] = newdist
    pred[v] = u
    end

    end
    end

    # retrieve path 
    v = argmax(x -> dist[x], vertices(g))
    push!(path, v)

    while pred[v] != v
    v = pred[v]
    push!(path, v)
    end

    # return 
    return reverse(path)
    =#

end

function getDijkstraRoute(data::SBRPData, consumed_times::Dict{Tuple{Int, Int}, Float64}, idxs_blocks::Vi)::Vi

    @debug "Get dijkstra route"

    # attrs

    @debug "Get params"

    m::Int = findlast(idx::Int -> any(i::Int -> consumed_times[(idx, i)] <= data.T, data.B[idx]), idxs_blocks)

    # edge case
    m == 0 && return tour

    tour::Vi       = Vi()
    B::VVi         = data.B
    A::Arcs        = data.D.A
    T::Float64     = data.T
    first_idx::Int = first(idxs_blocks)
    last_idx::Int  = idxs_blocks[m]

    # initialize tour
    
    @debug "Initialize tour"

    push!(tour, B[last_idx][findmin(i -> consumed_times[(last_idx, i)], B[last_idx])[2]])

    # BFS
    
    @debug "BFS"

    for i::Int in reverse(2:m)

        # get indexes

        @debug "Get curr and next block indexes"

        idx_curr_block::Int = idxs_blocks[i]
        idx_prev_block::Int = idxs_blocks[i - 1]

        # get curr and prev blocks

        @debug "Get curr and next block"

        curr_block::Vi = B[idx_curr_block]
        prev_block::Vi = B[idx_prev_block]

        # get intersection

        @debug "Get intersection"

        intersection::Vi = ∩(curr_block, prev_block)

        # get block time

        @debug "Get curr block time"

        curr_block_time::Float64 = blockTime(data, B[idx_curr_block])

        # if no intersecting nodes push backwards candidate

        @debug "No intersecting nodes edge-case"

        j::Int = last(tour)
        if !in(j, intersection)
            candidates::Vi = filter(
                                    i::Int -> consumed_times[(idx_prev_block, i)] + time(data, Arc(i, j)) + curr_block_time <= consumed_times[(idx_curr_block, j)], 
                                    prev_block
                                   )
            push!(tour, first(candidates))
        end

    end

    return reverse(tour)

end

"""
Countings
"""
N_FEASIBLE::Int        = 0
N_INFEASIBLE::Int      = 0
DECODING_TIME::Float64 = 0

function getIdxsBlocks(chromosome::Array{Float64}, data::SBRPData)::Vi
    # attrs
    m::Int = length(data.B)

    # get all the alleles with weight > 0.5 and sort them in reversed order 
    permutation::Array{Tuple{Float64, Int64}} = Array{Tuple{Float64, Int64}}(undef, m)

    for (idx::Int, key::Float64) in enumerate(chromosome)
        permutation[idx] = Tuple{Float64, Int}((key, idx))
    end

    sort!(permutation, rev = true)

    # consider only blocks with allele > 0.5
    #  filter!(s -> s[1] > 0.5, permutation)

    # return indexes of the selected blocks
    return map((key, idx)::Tuple{Float64, Int} -> idx, permutation)
end

function decode!(chromosome::Array{Float64}, data::SBRPData, rewrite::Bool)::Float64
    # import blogal variables
    global N_FEASIBLE
    global N_INFEASIBLE 
    global DECODING_TIME 

    DECODING_TIME += curr_decode_time::Float64 = @elapsed begin

        # attrs
        B::VVi = data.B

        # get min tour time
        idxs_blocks::Vi = getIdxsBlocks(chromosome, data)
        tour_time::Float64, consumed_times::Dict{Tuple{Int, Int}, Float64} = dijkstra(data, idxs_blocks)

        idxs_blocks = idxs_blocks[begin:findlast(idx::Int -> any(i::Int -> consumed_times[(idx, i)] <= data.T, data.B[idx]), idxs_blocks)]
    end

    if tour_time > data.T # check feasibility

        throw(InvalidStateException("Not here"))
        #      println("Infeasile with time ", tour_time)
        N_INFEASIBLE += 1
        return -∞

    else # return profit
        #      println("Feasile with time ", tour_time)
        N_FEASIBLE += 1
        return sum(idx::Int -> data.profits[B[idx]], idxs_blocks)

    end
end

function runBRKGAModel(data::SBRPData, app::Dict{String, Any})::Tuple{SBRPSolution, Dict{String, String}}

    global N_FEASIBLE
    global N_INFEASIBLE 
    global DECODING_TIME 

    N_FEASIBLE    = 0
    N_INFEASIBLE  = 0
    DECODING_TIME = 0

    # attrs
    B::VVi                     = data.B
    m::Int                     = length(B)
    verbose::Bool              = false
    info::Dict{String, String} = Dict{String, String}()

    ########################################
    #=
    profit, visited_blocks, tour = NearestNeighborhoodHeuristic.solve(data)
    visited_blocks_weights = Dict{Vi, Real}(visited_blocks[i] => 1 - i * 1e-3 for i in 1:length(visited_blocks))
    initial_chromosome = [block in keys(visited_blocks_weights) ? visited_blocks_weights[block] : 0.0 for block in B]

    tour = add_blocks(tour, visited_blocks)
    println(Data.SBRP.tour_time(data, tour, visited_blocks))

    # get blocks indexes
    idxs_blocks = get_idxs_blocks(initial_chromosome, data)
    println([findall(block′ -> block′ == block, B)[1] for block in visited_blocks])
    println(idxs_blocks)

    # get dijkstra time matrix
    tour_time, consumed_times = dijkstra(data, idxs_blocks)
    println(tour_time)

    for i in 2:length(idxs_blocks)
    idx_block, idx_prev_block = idxs_blocks[i], idxs_blocks[i - 1]
    curr_block_time = time_block(data, B[idx_block])
    println("ID block ", idx_block)
    for i in B[idx_block]
    rs = consumed_times[(idx_block, i)]
    println("\t", i, ": ", rs)
    for j in B[idx_prev_block]
    ls = consumed_times[(idx_prev_block, j)] + Data.SBRP.time(data, (j, i)) + curr_block_time
    if ls == rs
    println("\t\t", j, " is a predecessor, since ", ls, " == ", rs)
    else
    println("\t\t", j, " is not a predecessor, since ", ls, " != ", rs)
    end
    end
    end
    end

    # get tour
    tour = get_dijkstra_route(data, consumed_times, idxs_blocks)
    tour = add_blocks(tour, visited_blocks)
    println(Data.SBRP.tour_time(data, tour, visited_blocks))
    error("stop here")
    =#

    #=
    best_chromosome = [0.3589048242785371, 0.6319804041946304, 0.020200574214842337, 0.5084282197994525, 0.17132376723018838, 0.37440283851827494, 0.7616929330314515, 0.2051392112493251, 0.668209909295016, 0.1805155431288492, 0.605231415945501, 0.538706206753417, 0.643018709267497, 0.1052176490620893, 0.49865714857109733, 0.6203897012759658, 0.9649904947753347, 0.4654304539026679, 0.8277325292779287, 0.2817681302706594, 0.9441870264105543, 0.8417186847833293, 0.3046339513279499, 0.6076865663630291, 0.4873328865981925, 0.3277004411277693, 0.17926631433639972, 0.43547549746888614, 0.7786574742039711, 0.2726097587101328, 0.6273702912610448, 0.005491993320564159, 0.6610746153847575, 0.255518541943349, 0.4894043201742333, 0.11360851299666863, 0.7494590712211553, 0.881969426902371, 0.35435521386670077, 0.7235108062951003, 0.38944198819446996, 0.04768164110918938, 0.3324132325220943, 0.060793222727230756, 0.9712003842759238, 0.8285673944159864, 0.14717881496546692, 0.8982163374481262, 0.4816582024717766, 0.5132717889739846]  

    # get blocks indexes
    idxs_blocks = get_idxs_blocks(best_chromosome, data)

    # get dijkstra time matrix
    tour_time, consumed_times = dijkstra(data, idxs_blocks)

    for i in 2:length(idxs_blocks)
    idx_block, idx_prev_block = idxs_blocks[i], idxs_blocks[i - 1]
    curr_block_time = time_block(data, B[idx_block])
    println("ID block ", idx_block)
    for i in B[idx_block]
    rs = consumed_times[(idx_block, i)]
    println("\t", i, ": ", rs)
    for j in B[idx_prev_block]
    ls = consumed_times[(idx_prev_block, j)] + Data.SBRP.time(data, (j, i)) + curr_block_time
    if ls == rs
    println("\t\t", j, " is a predecessor, since ", ls, " == ", rs)
    else
    println("\t\t", j, " is not a predecessor, since ", ls, " != ", rs)
    end
    end
    end
    end

    # get tour
    tour = get_dijkstra_route(data, consumed_times, idxs_blocks)
    error("stop here")
    =#

    ########################################
    # Load configuration file and show basic info.
    ########################################

    conf_dir::String = string(app["brkga-conf"])
    conf::ConfParse   = ConfParse(conf_dir)
    parse_conf!(conf)

    seed::Int64                          = Base.parse(Int64, retrieve(conf, "seed"))
    stop_rule::StopRule                  = parseRule(StopRule, retrieve(conf, "stop_rule"))
    stop_argument::Union{Float64, Int64} = Base.parse(stop_rule == TARGET ? Float64 : Int64, retrieve(conf, "stop_argument"))
    maximum_time::Float64                = Base.parse(Float64, retrieve(conf, "maximum_time"))
    verbose                              = Base.parse(Bool, retrieve(conf, "verbose"))
    perform_evolution::Bool              = Base.parse(Bool, retrieve(conf, "evolution"))

    maximum_time <= 0.0 && ArgumentError("Maximum time must be larger than 0.0. Given $maximum_time.")

    ########################################
    # Load config file and show basic info.
    ########################################

    brkga_params, control_params = load_configuration(retrieve(conf, "brkg_params"))

    if verbose
        @info """
        ------------------------------------------------------
        > Experiment started at $(Dates.now())
        > Configuration: $conf_dir
        > Algorithm Parameters:
        """

        if !perform_evolution
            @info ">    - Simple multi-start: on (no evolutionary operators)"
        else
            # log BRKGA parameters
            for field in fieldnames(BrkgaParams)
                @info ">  - $field $(getfield(brkga_params, field))"
            end

            # log control parameters
            for field in fieldnames(ExternalControlParams)
                @info ">  - $field $(getfield(control_params, field))"
            end

            @info """
            > Seed: $seed
            > Stop rule: $stop_rule
            > Stop argument: $stop_argument
            > Maximum time (s): $maximum_time
            > Number of parallel threads for decoding: $(Threads.nthreads())
            ------------------------------------------------------"""
        end

        ########################################
        # Adjust BRKGA parameters
        ########################################

        @info "\n[$(Dates.Time(Dates.now()))] Generating initial tour..."
    end

    # Generate a greedy solution to be used as warm start for BRKGA.
    profit::Float64, visited_blocks::VVi, tour::Vi = solveNearestNeighborhood(data)
    verbose && @info "Initial profit: $profit"

    ########################################
    # Build the BRKGA data structures and initialize
    ########################################

    @info "\n[$(Dates.Time(Dates.now()))] Building BRKGA data..."

    # Usually, it is a good idea to set the population size
    # proportional to the instance size.
    brkga_params.population_size = min(brkga_params.population_size, 10 * m)
    @info "New population size: $(brkga_params.population_size)"

    # Chromosome size is the number of nodes.
    # Each chromosome represents a permutation of nodes.
    brkga_data = build_brkga(data, decode!, MAXIMIZE, seed, m, brkga_params, perform_evolution)

    # To inject the initial tour, we need to create chromosome representing that
    # solution. First, we create a set of keys to be used in the chromosome.
    # Then, we visit each required and profitable arc in the routes and assign to they keys.
    visited_blocks_weights::Dict{Vi, Float64} = Dict{Vi, Float64}(map(idx::Int -> visited_blocks[idx] => 1 - idx * 1e-3, 1:length(visited_blocks)))

    initial_chromosome::Vector{Float64} = map(block::Vi -> block in keys(visited_blocks_weights) ? visited_blocks_weights[block] : 0.0, B)

    # Inject the warm start solution in the initial population.
    set_initial_population!(brkga_data, [initial_chromosome])

    # NOTE: don't forget to initialize the algorithm.
    verbose && @info "\n[$(Dates.Time(Dates.now()))] Initializing BRKGA data..."
    initialize!(brkga_data)

    ########################################
    # Warm up the script/code
    ########################################

    # To make sure we are timing the runs correctly, we run some warmup
    # iterations with bogus data. Warmup is always recommended for script
    # languages. Here, we call the most used methods.
    verbose && @info "\n[$(Dates.Time(Dates.now()))] Warming up..."

    bogus_data = deepcopy(brkga_data)

    evolve!(bogus_data, 2)

    path_relink!(bogus_data, brkga_params.pr_type, brkga_params.pr_selection, (x, y) -> 1.0, (x, y) -> true, 0, 0.5, 1, 10.0, 1.0)

    best_cost::Float64 = get_best_fitness(brkga_data)

    best_chromosome::Vector{Float64} = get_best_chromosome(brkga_data)

    bogus_data = nothing

    ########################################
    # Evolving
    ########################################

    verbose && @info "\n[$(Dates.Time(Dates.now()))] Evolving..."
    verbose && @info "* Iteration | Cost | CurrentTime"

    iteration::Int              = 0
    last_update_time::Float64   = 0.0
    last_update_iteration::Int  = 0
    large_offset::Int           = 0
    path_relink_time::Float64   = 0.0
    num_path_relink_calls::Int  = 0
    num_homogenities::Int       = 0
    num_best_improvements::Int  = 0
    num_elite_improvements::Int = 0
    run::Bool                   = true
    start_time::Float64         = Base.time()

    # Main optimization loop. We evolve one generation at time,
    # keeping track of all changes during such process.
    while run
        iteration += 1

        # Evolves one iteration.
        evolve!(brkga_data)

        # Checks the current results and holds the best.
        fitness::Float64 = get_best_fitness(brkga_data)

        if fitness > best_cost

            last_update_time = Base.time() - start_time

            update_offset::Int = iteration - last_update_iteration

            if large_offset < update_offset 
                large_offset = update_offset
            end

            last_update_iteration = iteration
            best_cost             = fitness
            best_chromosome       = get_best_chromosome(brkga_data)

            verbose && @info @sprintf("* %d | %.0f | %.2f \n", iteration, best_cost, last_update_time)

        end

        iter_without_improvement::Int = iteration - last_update_iteration

        # Here, we call the path relink when the algorithm gets stuck for
        # `exchange_interval` iterations. Obviously, we can use many other ways
        # of hybridization.
        if control_params.exchange_interval > 0 &&
            iter_without_improvement > 0 &&
            (iter_without_improvement % control_params.exchange_interval == 0)

            verbose && @info "Performing path relink at $iteration..."
            num_path_relink_calls += 1

            pr_now::Float64 = Base.time()
            result = path_relink!(
                                  brkga_data,
                                  brkga_params.pr_type,
                                  brkga_params.pr_selection,
                                  kendall_tau_distance,
                                  affect_solution_kendall_tau,
                                  brkga_params.pr_number_pairs,
                                  brkga_params.pr_minimum_distance,
                                  1, # block_size doesn't matter for permutation.
                                  maximum_time - (Base.time() - start_time),
                                  brkga_params.pr_percentage
                                 )

            pr_time::Float64 = Base.time() - pr_now
            path_relink_time += pr_time

            if result == TOO_HOMOGENEOUS
                num_homogenities += 1
                verbose && @info "- Populations are too too homogeneous | " *
                "Elapsed time: $(@sprintf("%.2f", pr_time))"

            elseif result == NO_IMPROVEMENT
                verbose && @info "- No improvement found | " *
                "Elapsed time: $(@sprintf("%.2f", pr_time))"

            elseif result == ELITE_IMPROVEMENT
                num_elite_improvements += 1
                verbose && @info "- Improvement on the elite set but " *
                "not in the best individual | " *
                "Elapsed time: $(@sprintf("%.2f", pr_time))"

            elseif result == BEST_IMPROVEMENT

                num_best_improvements += 1

                fitness = get_best_fitness(brkga_data)

                verbose && @info "- Best individual improvement: $fitness | " *
                "Elapsed time: $(@sprintf("%.2f", pr_time))"

                if fitness < best_cost

                    last_update_time      = Base.time() - start_time
                    update_offset         = iteration - last_update_iteration

                    if large_offset < update_offset 
                        large_offset = update_offset
                    end

                    last_update_iteration = iteration
                    best_cost             = fitness
                    best_chromosome       = get_best_chromosome(brkga_data)

                    verbose && @info @sprintf("* %d | %.0f | %.2f \n", iteration, best_cost, last_update_time)
                end
            end
        end

        # Check stop criteria.
        run = !(
                Base.time() - start_time > maximum_time ||
                (stop_rule == GENERATIONS && Float64(iteration) == stop_argument) ||
                (stop_rule == IMPROVEMENT && Float64(iter_without_improvement) >= stop_argument) ||
                (stop_rule == TARGET && best_cost <= stop_argument)
               )
    end

    total_elapsed_time::Float64 = Base.time() - start_time
    total_num_iterations::Int = iteration

    if verbose
        @info "[$(Dates.Time(Dates.now()))] End of optimization\n"
        @info "Total number of iterations: $total_num_iterations\n"
        @info "Last update iteration: $last_update_iteration\n"
        @info @sprintf("Total optimization time: %.2f\n", total_elapsed_time)
        @info @sprintf("Last update time: %.2f\n", last_update_time)
        @info "Large number of iterations between improvements: $large_offset\n"
        @info @sprintf("Total path relink time: %.2f\n", path_relink_time)
        @info "Total path relink calls: $num_path_relink_calls\n"
        @info "Number of homogenities: $num_homogenities\n"
        @info "Improvements in the elite set: $num_elite_improvements\n"
        @info "Best individual improvements: $num_best_improvements\n"
    end

    ########################################
    # Extracting the solution
    ########################################

    # get blocks indexes
    idxs_blocks::Vi = getIdxsBlocks(best_chromosome, data)

    # get dijkstra time matrix
    tour_time::Float64, consumed_times::Dict{Tuple{Int, Int}, Float64} = dijkstra(data, idxs_blocks)

    # get tour
    tour = getDijkstraRoute(data, consumed_times, idxs_blocks)

    @info @sprintf("Last index %d", findlast(idx::Int -> any(i::Int -> consumed_times[(idx, i)] <= data.T, data.B[idx]), idxs_blocks))

    # update blocks indexes
    idxs_blocks = idxs_blocks[1:findlast(idx::Int -> any(i::Int -> consumed_times[(idx, i)] <= data.T, data.B[idx]), idxs_blocks)]

    # get visited blocks
    blocks::VVi = map(idx::Int -> B[idx], idxs_blocks)

    # log
    #  info["cost"], info["solverTime"] = string(∑(data.profits[block] for block in blocks)), string(total_elapsed_time)
    info["cost"]       = @sprintf("%.2f", best_cost)
    info["solverTime"] = string(total_elapsed_time)
    info["numVisitedBlocks"] = string(length(blocks))

    @info @sprintf("BRKGA cost: %.2f", best_cost)
    @info @sprintf("Blocks cost: %.2f", sum(block::Vi -> data.profits[block], blocks))
    @info @sprintf("Calculated cost: %s", info["cost"])
    
    @info "Return"

    return SBRPSolution(tour, blocks), info

end
