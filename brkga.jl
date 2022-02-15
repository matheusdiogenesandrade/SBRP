module BRKGA

include("symbols.jl")

export run_brkga

using ..Data
using ..Data.SBRP
using ..NearestNeighborhoodHeuristic
using ..Solution

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
function parse(::Type{StopRule}, value::String)::StopRule
    local_value = uppercase(strip(value)[1])
    dict = Dict{Char, StopRule}('G' => GENERATIONS, 'T' => TARGET, 'I' => IMPROVEMENT)

    local_value in keys(dict) && return dict[local_value]
    throw(ArgumentError("cannot parse $value as StopRule"))
end

function dijkstra(data::SBRPData, idxs_blocks::Vi)

  m = length(idxs_blocks)

  # edge case
  m == 0 && return 0.0, nothing

  # attrs
  B, A, T, first_idx = data.B, data.D.A, data.T, first(idxs_blocks)

  # dijkstra distance <<idx_block, node>, min path time> stores the min path time from the depot to the node `node' in the block at `idx_block'
  consumed_times = Dict{Tuple{Int, Int}, Real}((idx_block, i) => T + 1 for idx_block in idxs_blocks for i in B[idx_block])

  # initial times
  first_block_time = time_block(data, B[first_idx])
  [consumed_times[(first_idx, i)] = first_block_time for i in B[first_idx]]

  # BFS
  for i in 1:m - 1
    # get indexes
    idx_curr_block, idx_next_block = idxs_blocks[i], idxs_blocks[i + 1]

    # get curr and next blocks
    curr_block, next_block = B[idx_curr_block], B[idx_next_block]

    # calculate nest block service time
    next_block_time = time_block(data, next_block)

    # intersecting nodes
    intersection = Si(∩(curr_block, next_block))
    for j in intersection
      consumed_times[(idx_next_block, j)] = next_block_time + consumed_times[(idx_curr_block, j)]
    end

    # non intersecting nodes
    for i in curr_block
      for j in next_block
        (i, j) in keys(data.D.distance) && (consumed_times[(idx_next_block, j)] = min(consumed_times[(idx_next_block, j)], consumed_times[(idx_curr_block, i)] + Data.SBRP.time(data, (i, j)) + next_block_time))
      end
    end
  end
  
  return min([consumed_times[(last(idxs_blocks), i)] for i in B[last(idxs_blocks)]]...), consumed_times
end

function get_dijkstra_route(data::SBRPData, consumed_times::Dict{Tuple{Int, Int}, Real}, idxs_blocks::Vi)

  tour, m = Vi(), length(idxs_blocks)

  # edge case
  m == 0 && return tour

  # attrs
  B, A, T, first_idx, last_idx = data.B, data.D.A, data.T, first(idxs_blocks), last(idxs_blocks)

  # initialize tour
  
  push!(tour, findmin([(consumed_times[(last_idx, i)], i) for i in B[last_idx]])[1][2])

  # BFS
  for i in reverse(2:m)
    # get indexes
    idx_curr_block, idx_prev_block = idxs_blocks[i], idxs_blocks[i - 1]

    # get curr and prev blocks
    curr_block, prev_block = B[idx_curr_block], B[idx_prev_block]

    # get intersection
    intersection = ∩(curr_block, prev_block)

    # get block time
    curr_block_time = time_block(data, B[idx_curr_block])

#    println("In block ", idx_curr_block, " previous node ", last(tour))

    # if no intersecting nodes push backwards candidate
    !in(last(tour), intersection) && push!(tour, first([i for (i, j) in δ⁻(A, last(tour)) if i in prev_block && consumed_times[(idx_prev_block, i)] + Data.SBRP.time(data, (i, j)) + curr_block_time <= consumed_times[(idx_curr_block, j)]]))

  end

  return reverse(tour)
   
end

"""
Countings
"""
N_FEASIBLE, N_INFEASIBLE, DECODING_TIME = 0, 0, 0.0

function get_idxs_blocks(chromosome::Array{Float64}, data::SBRPData)
  # attrs
  m = length(data.B)

  # get all the alleles with weight > 0.5 and sort them in reversed order 
  permutation = Array{Tuple{Float64, Int64}}(undef, m)
  [permutation[idx] = (key, idx) for (idx, key) in enumerate(chromosome)]
  sort!(permutation, rev = true)

  # consider only blocks with allele > 0.5
  filter!(s -> s[1] > 0.5, permutation)

  # return indexes of the selected blocks
  return [idx for (key, idx) in permutation]
end

function decode!(chromosome::Array{Float64}, data::SBRPData, rewrite::Bool)::Float64
  # inport blogal variables
  global N_FEASIBLE
  global N_INFEASIBLE 
  global DECODING_TIME 

  DECODING_TIME += @elapsed begin

    # attrs
    B = data.B

    # get min tour time
    idxs_blocks = get_idxs_blocks(chromosome, data)
    tour_time, consumed_times = dijkstra(data, idxs_blocks)

    
    if tour_time > data.T # check feasibility

#      println("Infeasile with time ", tour_time)
      N_INFEASIBLE += 1
      return -∞

    else # return profit

#      println("Feasile with time ", tour_time)
      N_FEASIBLE += 1
      return ∑(data.profits[B[idx]] for idx in idxs_blocks)

    end
  end
end

function run_brkga(conf_dir::String, data::SBRPData)

  global N_FEASIBLE
  global N_INFEASIBLE 
  global DECODING_TIME 

  N_FEASIBLE, N_INFEASIBLE, DECODING_TIME = 0, 0, 0.0

  # attrs
  B = data.B
  m = length(B)
  verbose, info = false, Dict{String, String}()

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

  conf              = ConfParse(conf_dir)
  parse_conf!(conf)

  seed              = Base.parse(Int64, retrieve(conf, "seed"))
  stop_rule         = parse(StopRule, retrieve(conf, "stop_rule"))
  stop_argument     = Base.parse(stop_rule == TARGET ? Float64 : Int64, retrieve(conf, "stop_argument"))
  maximum_time      = Base.parse(Float64, retrieve(conf, "maximum_time"))
  verbose           = Base.parse(Bool, retrieve(conf, "verbose"))
  perform_evolution = Base.parse(Bool, retrieve(conf, "evolution"))

  maximum_time <= 0.0 && error("Maximum time must be larger than 0.0. Given $maximum_time.")

  ########################################
  # Load config file and show basic info.
  ########################################

  brkga_params, control_params = load_configuration(retrieve(conf, "brkg_params"))

  verbose && flush_println("""
    ------------------------------------------------------
    > Experiment started at $(Dates.now())
    > Configuration: $conf_dir
    > Algorithm Parameters:
    """)

  if !perform_evolution
    verbose && println(">    - Simple multi-start: on (no evolutionary operators)")
  else
    # log BRKGA parameters
    [(verbose && flush_println(">  - $field $(getfield(brkga_params, field))")) for field in fieldnames(BrkgaParams)]
    # log control parameters
    [(verbose && flush_println(">  - $field $(getfield(control_params, field))")) for field in fieldnames(ExternalControlParams)]
    verbose && flush_println("""
            > Seed: $seed
            > Stop rule: $stop_rule
            > Stop argument: $stop_argument
            > Maximum time (s): $maximum_time
            > Number of parallel threads for decoding: $(Threads.nthreads())
            ------------------------------------------------------""")
  end

  ########################################
  # Adjust BRKGA parameters
  ########################################

  verbose && flush_println("\n[$(Dates.Time(Dates.now()))] Generating initial tour...")

  # Generate a greedy solution to be used as warm start for BRKGA.
  profit, visited_blocks, tour = NearestNeighborhoodHeuristic.solve(data)
  verbose && flush_println("Initial profit: $profit")

  ########################################
  # Build the BRKGA data structures and initialize
  ########################################

  flush_println("\n[$(Dates.Time(Dates.now()))] Building BRKGA data...")

  # Usually, it is a good idea to set the population size
  # proportional to the instance size.
  brkga_params.population_size = min(brkga_params.population_size, 10 * m)
  flush_println("New population size: $(brkga_params.population_size)")

  # Chromosome size is the number of nodes.
  # Each chromosome represents a permutation of nodes.
  brkga_data = build_brkga(data, decode!, MAXIMIZE, seed, m, brkga_params, perform_evolution)

  # To inject the initial tour, we need to create chromosome representing that
  # solution. First, we create a set of keys to be used in the chromosome.
  # Then, we visit each required and profitable arc in the routes and assign to they keys.
  visited_blocks_weights = Dict{Vi, Real}(visited_blocks[i] => 1 - i * 1e-3 for i in 1:length(visited_blocks))
  initial_chromosome = [block in keys(visited_blocks_weights) ? visited_blocks_weights[block] : 0.0 for block in B]
#  initial_chromosome = [(rand()%0.5) + (block in visited_blocks ?  0.5 : 0.0) for block in B]


  # Inject the warm start solution in the initial population.
  set_initial_population!(brkga_data, [initial_chromosome])

  # NOTE: don't forget to initialize the algorithm.
  verbose && flush_println("\n[$(Dates.Time(Dates.now()))] Initializing BRKGA data...")
  initialize!(brkga_data)

  ########################################
  # Warm up the script/code
  ########################################

  # To make sure we are timing the runs correctly, we run some warmup
  # iterations with bogus data. Warmup is always recommended for script
  # languages. Here, we call the most used methods.
  verbose && flush_println("\n[$(Dates.Time(Dates.now()))] Warming up...")

  bogus_data = deepcopy(brkga_data)
  evolve!(bogus_data, 2)
  path_relink!(bogus_data, brkga_params.pr_type, brkga_params.pr_selection, (x, y) -> 1.0, (x, y) -> true, 0, 0.5, 1, 10.0, 1.0)
  get_best_fitness(brkga_data)
  get_best_chromosome(brkga_data)
  bogus_data = nothing

  ########################################
  # Evolving
  ########################################

  verbose && flush_println("\n[$(Dates.Time(Dates.now()))] Evolving...")
  verbose && flush_println("* Iteration | Cost | CurrentTime")

#  best_cost = -Inf
  best_cost = profit
  best_chromosome = initial_chromosome

  iteration = 0
  last_update_time = 0.0
  last_update_iteration = 0
  large_offset = 0
  path_relink_time = 0.0
  num_path_relink_calls = 0
  num_homogenities = 0
  num_best_improvements = 0
  num_elite_improvements = 0
  run = true
  start_time = time()

  # Main optimization loop. We evolve one generation at time,
  # keeping track of all changes during such process.
  while run
    iteration += 1

    # Evolves one iteration.
    evolve!(brkga_data)

    # Checks the current results and holds the best.
    fitness = get_best_fitness(brkga_data)
    if fitness > best_cost
      last_update_time = time() - start_time
      update_offset = iteration - last_update_iteration

      large_offset < update_offset && (large_offset = update_offset)

      last_update_iteration = iteration
      best_cost = fitness
      best_chromosome = get_best_chromosome(brkga_data)

      verbose && @printf("* %d | %.0f | %.2f \n", iteration, best_cost, last_update_time)
      
    end

    iter_without_improvement = iteration - last_update_iteration

    # Here, we call the path relink when the algorithm gets stuck for
    # `exchange_interval` iterations. Obviously, we can use many other ways
    # of hybridization.
    if control_params.exchange_interval > 0 &&
      iter_without_improvement > 0 &&
      (iter_without_improvement % control_params.exchange_interval == 0)

      verbose && flush_println("Performing path relink at $iteration...")
      num_path_relink_calls += 1

      pr_now = time()
      result = path_relink!(
                            brkga_data,
                            brkga_params.pr_type,
                            brkga_params.pr_selection,
                            kendall_tau_distance,
                            affect_solution_kendall_tau,
                            brkga_params.pr_number_pairs,
                            brkga_params.pr_minimum_distance,
                            1, # block_size doesn't matter for permutation.
                            maximum_time - (time() - start_time),
                            brkga_params.pr_percentage
                           )

      pr_time = time() - pr_now
      path_relink_time += pr_time

      if result == TOO_HOMOGENEOUS
        num_homogenities += 1
        verbose && flush_println("- Populations are too too homogeneous | " *
                "Elapsed time: $(@sprintf("%.2f", pr_time))")

      elseif result == NO_IMPROVEMENT
        verbose && flush_println("- No improvement found | " *
                "Elapsed time: $(@sprintf("%.2f", pr_time))")

      elseif result == ELITE_IMPROVEMENT
        num_elite_improvements += 1
        verbose && flush_println("- Improvement on the elite set but " *
                "not in the best individual | " *
                "Elapsed time: $(@sprintf("%.2f", pr_time))")

      elseif result == BEST_IMPROVEMENT
        num_best_improvements += 1
        fitness = get_best_fitness(brkga_data)
        verbose && flush_println("- Best individual improvement: $fitness | " *
                "Elapsed time: $(@sprintf("%.2f", pr_time))")
        if fitness < best_cost
          last_update_time = time() - start_time
          update_offset = iteration - last_update_iteration

          large_offset < update_offset && (large_offset = update_offset)

          last_update_iteration = iteration
          best_cost = fitness
          best_chromosome = get_best_chromosome(brkga_data)

          verbose && @printf("* %d | %.0f | %.2f \n", iteration, best_cost,
                  last_update_time)
        end
      end
    end

    # Check stop criteria.
    run = !(
            (time() - start_time > maximum_time) ||
            (stop_rule == GENERATIONS && Float64(iteration) == stop_argument) ||
            (stop_rule == IMPROVEMENT &&
             Float64(iter_without_improvement) >= stop_argument) ||
            (stop_rule == TARGET && best_cost <= stop_argument)
           )
  end
  total_elapsed_time = time() - start_time
  total_num_iterations = iteration

  if verbose
    flush_println("[$(Dates.Time(Dates.now()))] End of optimization")
    flush_println()
    flush_println("Total number of iterations: $total_num_iterations")
    flush_println("Last update iteration: $last_update_iteration")
    @printf("Total optimization time: %.2f\n", total_elapsed_time)
    @printf("Last update time: %.2f\n", last_update_time)
    flush_println("Large number of iterations between improvements: $large_offset")

    @printf("Total path relink time: %.2f\n", path_relink_time)
    flush_println("Total path relink calls: $num_path_relink_calls")
    flush_println("Number of homogenities: $num_homogenities")
    flush_println("Improvements in the elite set: $num_elite_improvements")
    flush_println("Best individual improvements: $num_best_improvements")
  end

  ########################################
  # Extracting the solution
  ########################################

  # get blocks indexes
  idxs_blocks = get_idxs_blocks(best_chromosome, data)

  # get dijkstra time matrix
  tour_time, consumed_times = dijkstra(data, idxs_blocks)

  # get tour
  tour = get_dijkstra_route(data, consumed_times, idxs_blocks)

  # get visited blocks
  blocks = [B[idx_block] for idx_block in idxs_blocks]

  # log
  info["cost"], info["solverTime"] = string(∑(data.profits[block] for block in blocks)), string(total_elapsed_time)

  return tour, info, blocks

end

end
