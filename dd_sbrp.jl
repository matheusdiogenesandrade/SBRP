module DD_SBRP

include("symbols.jl")
include("dd.jl")

export run_dd

using ..Data
using ..Data.SBRP
using ..NearestNeighborhoodHeuristic
using ..Solution
using .DecisionDiagram

using Dates

function run_dd(data::SBRPData)
  # setup
  T, V, depot, B = data.T, data.D.V, data.depot, data.B 
  info = Dict{String, String}()

  # helpers
  node_blocks = Dict{Int, VVi}(i => filter(block -> i in block, B) for i in keys(V))
  BlockDistance = Dict{Int, Float64}
  get_intersections(nodes::Si) = [Si()]

  #=
  Building DD
  =#
  # parameters
  variables = Vector{Any}(1:length(B))
  domain = Vector{Any}(variables)
  initial_state = State("times" => BlockDistance(depot => 0.0), "invalid_nodes" => Si([depot]), "remaining_blocks" => Si(variables), "profit" => 0.0, "serviced_block" => nothing) 
  empty_domain = -1 

  # callbacks
  function get_next_state(idx_variable::Int, state::DecisionDiagram.State, candidate::Any, dd::DecisionDiagram.Instance)
    # edge case
    candidate == empty_domain && return [State()]

    # next block
    block = B[candidate]

    # times
    times′ = BlockDistance(i => min((time + Data.SBRP.time(data, (i, j)) for (j, time) in state["times"])...) + Data.SBRP.time_block(data, block) for i in block if !in(i, state["invalid_nodes"]))
    filter!(p -> p[2] <= T, times′)

    # edge case
    isempty(times′)  && return [State()]

    # get extra nodes to exclude
    excludeds = get_intersections(Si(keys(times′)))

    # states
    states = []

    for excluded in excludeds

      state′ = State()
      state′["invalid_nodes"] = union(state["invalid_nodes"], excluded)
      state′["times"] = filter(p -> !in(p[1], excluded), times′)

      # edge case
      isempty(state′["times"])  && return [State()]

      state′["remaining_blocks"] = setdiff(state["remaining_blocks"], candidate)
      state′["profit"] = state["profit"] + data.profits[block]
      state′["serviced_block"] = candidate

      length(state′["times"]) > 0 && push!(states, state′)
    end

    return states
  end

  get_candidates(idx_variable::Int, state::DecisionDiagram.State, dd::DecisionDiagram.Instance) = state["remaining_blocks"]

  function not_redundant_state(idx_variable::Int, state::State, states::Vector{State}) 
    isempty(state) && return false
    
    all([time > T for (i, time) in state["times"]]) && return true

    redundant = false

    for state′ in states
      # edge case
      state′ == state && continue

      # if not same remaining blocks 
      state′["remaining_blocks"] != state["remaining_blocks"] && continue

      # if state′ profit is worst than state
      state′["profit"] < state["profit"] && continue

      # get intersection 
      intersection = intersection(keys(state′["times"]), keys(state["times"]))

      # is no intersection just get out
      isempty(intersection) && continue

      # if state′ "covers" state
      redundant = all([state′["times"][i] <= state["times"] for i in intersection])

      redundant && break

    end

    return !redundant 
  end
#  not_redundant_state(idx_variable::Int, state::State, states::Vector{State}) = true

  dd = DecisionDiagram.Instance(empty_domain, domain, variables, initial_state, get_candidates, get_next_state, not_redundant_state)
#  dd.LOG = true
  total_elapsed_time = @elapsed DecisionDiagram.run(dd)

  # get best solution
  best_var, best_state = 0, initial_state
  for idx_var in 1:length(dd.variables)
    for state in dd.states[idx_var]
      if state["profit"] > best_state["profit"]
        best_var, best_state = idx_var, state
      end
    end
  end

#  println(best_var)
#  println(best_state)

  ########################################
  # Extracting the solution
  ########################################
  # helpers
  find_valid_prev_nodes(node::Int, time_node::Float64, state::State) = keys(filter(p -> p[2] + Data.SBRP.time(data, (p[1], node)) <= time_node, state["times"]))
  find_valid_prev_states(node::Int, time_node::Float64, states::Vector{State}) = filter(state -> !isempty(find_valid_prev_nodes(node, time_node, state)), states)

  # edge case
  isempty(find_valid_prev_nodes(depot, T, best_state)) && error("A best solution is not feasible")

  # setup
  state = best_state
  node = first(find_valid_prev_nodes(depot, T, state))
  time = state["times"][node] - Data.SBRP.time_block(data, B[state["serviced_block"]])
  tour = [node, depot]

#  println("Tail state:\n\t$state")

  # traverse variables
  for idx_var in reverse(2:best_var)  
    # get valid previous states
    states = find_valid_prev_states(node, time, dd.states[idx_var - 1])

    # edge case
    isempty(states) && error("At variable $idx_var there is no valid previous state")

    # get previous state
    state = first(states)

#    println("$idx_var-th state:\n\t$state")

    # get previous nodes
    prev_nodes = find_valid_prev_nodes(node, time, state)

    # edge case
    isempty(prev_nodes) && error("There is no valid previous nodes")

    # get previous node
    prev_node = first(prev_nodes)

    # update 
    if prev_node != node
      node = prev_node
      prepend!(tour, node)
    end
    time = state["times"][node] - Data.SBRP.time_block(data, B[state["serviced_block"]])
  end

  for (idx_var, states) in enumerate(dd.states)
    println("Depth - $idx_var with $(length(states)) states") 
  end

  # add depot
  prepend!(tour, depot)
#  println(tour)

  # get serviced blocks
  blocks = [B[idx_block] for idx_block in setdiff(variables, best_state["remaining_blocks"])]

  # get total profit
#  println(sum(data.profits[block] for block in blocks))
  ∑(data.profits[block] for block in blocks) != best_state["profit"] && error("The profit does not match $(∑(data.profits[block] for block in blocks)) != $(best_state["profit"])")

  # log
  info["cost"], info["solverTime"] = string(best_state["profit"]), string(total_elapsed_time)

  return tour, info, blocks

end

end
