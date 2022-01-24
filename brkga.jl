module BRKGA

include("symbols.jl")

export run_brkga

using ..Data
using ..Data.SBRP
using ..NearestNeighborhoodHeuristic

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
  m == 0 && return 0.0

  # attrs
  B, A, T, first_idx = data.B, data.D.A, data.T, first(idxs_blocks)

  # dijkstra distance <<idx_block, node>, min path time> stores the min path time from the depot to the node `node' in the block at `idx_block'
  consumed_times = Dict{Tuple{Int, Int}, Real}((idx_block, i) => T + 1 for idx_block in idxs_blocks for i in B[idx_block])

  # initial times
  first_block_time = time_block(data, B[first_idx])
  [consumed_times[(first_idx, i)] = first_block_time for i in B[first_idx]]

  # BFS
  for i in 1:m - 1
    idx_curr_block = idxs_blocks[i]
    idx_next_block = idxs_blocks[i + 1]

    # get curr and next blocks
    curr_block, next_block = B[idx_curr_block], B[idx_next_block]

    # intersecting nodes
    for i in ∩(curr_block, next_block)
      consumed_times[(idx_next_block, i)] = min(consumed_times[(idx_next_block, i)], consumed_times[(idx_curr_block, i)])
    end

    # non intersecting nodes
    for (i, j) in δ⁺(A, curr_block)
      j in next_block && (consumed_times[(idx_next_block, j)] = min(consumed_times[(idx_next_block, j)], consumed_times[(idx_curr_block, i)] + Data.SBRP.time(data, (i, j))))
    end
  end
  
  return min([consumed_times[(last(idxs_blocks), i)] for i in B[last(idxs_blocks)]]...)
end

"""
Countings
"""
N_FEASIBLE = 0
N_INFEASIBLE = 0
DECODING_TIME = 0.0

function decode!(chromosome::Array{Float64}, data::SBRPData, rewrite::Bool)::Float64
  # inport blogal variables
  global N_FEASIBLE
  global N_INFEASIBLE 
  global DECODING_TIME 

  DECODING_TIME += @elapsed begin

    # attrs
    B = data.B
    m, tour_time = length(B), 0.0

    # get all the alleles with weight > 0.5 and sort them in reversed order 
    permutation = Array{Tuple{Float64, Int64}}(undef, m)
    [permutation[idx] = (key, idx) for (idx, key) in enumerate(chromosome)]
    sort!(permutation, rev = true)

    # consider only blocks with allele > 0.5
    filter!(s -> s[1] > 0.5, permutation)

    # get min tour time
    idxs_blocks = [idx for (key, idx) in permutation]
    tour_time = dijkstra(data, idxs_blocks)

    #  flush_println(idxs_blocks, " ", tour_time, " ", sum(data.profits[B[idx]] for idx in idxs_blocks))

    # check feasibility
    if tour_time > data.T

      N_INFEASIBLE += 1
      return -∞

    else # return profit

      N_FEASIBLE += 1
      return ∑(data.profits[B[idx]] for idx in idxs_blocks)

    end
  end
end

function run_brkga(app::Dict{String, Any}, data::SBRPData)
  #= 
  # params
  seed, configuration_file, num_generations = parse(Int64, app["brkga-seed"]), app["brkga-conf"], parse(Int, app["krkga-ngenerations"])

  seed              = parse(Int64, retrieve(conf, "seed"))

  # build 
  brkga_data, control_params = build_brkga(data, decode!, MAXIMIZE, seed, length(data.B), configuration_file)

  # init
  initialize!(brkga_data)

  # ...
  evolve!(brkga_data, num_generations)

  best_cost = get_best_fitness(brkga_data)

  @show best_cost
  =#
  B = data.B
  m = length(B)

  verbose = false 

  ########################################
  # Load configuration file and show basic info.
  ########################################

  conf              = ConfParse(app["brkga-conf"])
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
    > Configuration: $(app["brkga-conf"])
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
  initial_chromosome = [(rand()%0.5) + (block in visited_blocks ?  0.5 : 0.0) for block in B]

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

  best_cost = -Inf
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
    if fitness < best_cost
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
  
#  "$(@sprintf("%.2f", total_elapsed_time)) & " *
end

end
