module ExactDecisionDiagramTests

include("dd.jl")

using .DecisionDiagram
using Test

function knapsack_test(I, L)
  weights = [I[i][1] for i in 1:length(I)]
  #=
  Building DD
  =#
  # parameters
  empty_domain = -1                            # any value is valid, since the domain 0 always will be available (check get_candidates function definition), this happens because the problem modeling guarantees this.
  domain = Array{Any}([0, 1])                  # 0, if we do not select item, and 1 otherwise
  variables = Array{Any}(weights)              # the items weights
  initial_state = State("total_weight" => 0.0) # empty knapsack

  # callbacks
  get_next_state(idx_variable::Int64, state::DecisionDiagram.State, candidate::Any, dd::DecisionDiagram.Instance) = State("total_weight" => state["total_weight"] + (candidate == 1 ? dd.variables[idx_variable] : 0.0)) 

  get_candidates(idx_variable::Int64, state::DecisionDiagram.State, dd::DecisionDiagram.Instance) = state["total_weight"] + dd.variables[idx_variable] <= L ? [0, 1] : [0] 

  not_redundant_state(idx_variable::Int64, state::State, states::Array{State}) = true

  dd = DecisionDiagram.Instance(empty_domain, domain, variables, initial_state, get_candidates, get_next_state, not_redundant_state)
  # dd.LOG = true
  DecisionDiagram.run(dd)

  return dd
end

function knapsack_test_case_1()
  I, L = [(5, 3), (5, 5), (2, 2), (7, 4), (1, 1)], 7
  dd = knapsack_test(I, L)
  required_states = [
                     [0.0, 5.0], 
                     [0.0, 5.0], 
                     [0.0, 2.0, 5.0, 7.0], 
                     [0.0, 7.0, 2.0, 5.0], 
                     [0.0, 1.0, 7.0, 2.0, 3.0, 5.0, 6.0]
                    ]
  for idx in 1:length(dd.variables)
    states = [state["total_weight"] for state in dd.states[idx]]
    @test length(required_states[idx]) == length(states) && all([i in states for i in required_states[idx]])
  end
end

function knapsack_test_case_2()
  I, L = [(5, 3), (5, 5), (2, 2), (7, 4), (1, 1)], 10
  dd = knapsack_test(I, L)
  required_states = [
                     [0.0, 5.0], 
                     [0.0, 5.0, 10.0], 
                     [0.0, 2.0, 5.0, 7.0, 10.0], 
                     [0.0, 2.0, 9.0, 7.0, 5.0, 10.0], 
                     [0.0, 1.0, 9.0, 7.0, 2.0, 3.0, 5.0, 6.0, 10.0, 8.0]
                    ]
  for idx in 1:length(dd.variables)
    states = [state["total_weight"] for state in dd.states[idx]]
    @test length(required_states[idx]) == length(states) && all([i in states for i in required_states[idx]])
  end
end

function set_covering_test(E, S)
  #=
  Building DD
  =#
  # parameters
  empty_domain = -1                            # any value is valid, since the domain 0 always will be available (check get_candidates function definition), this happens because the problem modeling guarantees this.
  domain = Array{Any}([0, 1])                  # 0, if we do not select item, and 1 otherwise
  variables = Array{Any}(E)                    # the elements
  initial_state = State("uncovered_sets" => S) # all sets are uncovered

  # callbacks
  get_next_state(idx_variable::Int64, state::DecisionDiagram.State, candidate::Any, dd::DecisionDiagram.Instance) = candidate == 1 ? State("uncovered_sets" => [s for s in state["uncovered_sets"] if !in(dd.variables[idx_variable], s)]) : state

  get_candidates(idx_variable::Int64, state::DecisionDiagram.State, dd::DecisionDiagram.Instance) = isempty(findall(s -> last(s) == dd.variables[idx_variable], state["uncovered_sets"])) ? [0, 1] : [1] 

  not_redundant_state(idx_variable::Int64, state::State, states::Array{State}) = true

  dd = DecisionDiagram.Instance(empty_domain, domain, variables, initial_state, get_candidates, get_next_state, not_redundant_state)
  # dd.LOG = true
  DecisionDiagram.run(dd)

  return dd
end

function set_covering_test_case_1()
  E, S = 1:6, [[1, 2, 3], [1, 4, 5], [2, 4, 6]] # for this model, it is necessary that E and all elements of S are sorted.
  
  dd = set_covering_test(E, S) 
  first, second, third = S[1], S[2], S[3]
  required_states = [
                     [[first, second, third], [third]],
                     [[first, second, third], [second], [third], []],
                     [[second, third], [second], [third], []],
                     [[second, third], [second], [third], []],
                     [[third], []],
                     [[]],
                    ]
  for idx in 1:length(dd.variables)
    states = [state["uncovered_sets"] for state in dd.states[idx]]
    @test length(required_states[idx]) == length(states) && all([i in states for i in required_states[idx]])
  end
end

function packing_set_test(E, S)
  #=
  Building DD
  =#
  # parameters
  empty_domain = -1                            # any value is valid, since the domain 0 always will be available (check get_candidates function definition), this happens because the problem modeling guarantees this.
  domain = Array{Any}([0, 1])                  # 0, if we do not select item, and 1 otherwise
  variables = Array{Any}(E)                    # the elements
  initial_state = State("uncovered_sets" => S) # all sets are uncovered

  # callbacks
  get_next_state(idx_variable::Int64, state::DecisionDiagram.State, candidate::Any, dd::DecisionDiagram.Instance) = State("uncovered_sets" => candidate == 1 ? [s for s in state["uncovered_sets"] if !in(dd.variables[idx_variable], s)] : filter(s -> last(s) != dd.variables[idx_variable], state["uncovered_sets"]))

  get_candidates(idx_variable::Int64, state::DecisionDiagram.State, dd::DecisionDiagram.Instance) = isempty(findall(s -> !in(s, state["uncovered_sets"]) && dd.variables[idx_variable] in s, dd.initial_state["uncovered_sets"])) ? [0, 1] : [0] 

  not_redundant_state(idx_variable::Int64, state::State, states::Array{State}) = true

  dd = DecisionDiagram.Instance(empty_domain, domain, variables, initial_state, get_candidates, get_next_state, not_redundant_state)
  # dd.LOG = true
  DecisionDiagram.run(dd)

  return dd
end

function packing_set_test_case_1()
  E, S = 1:6, [[1, 2, 3], [1, 4, 5], [2, 4, 6]] # for this model, it is necessary that E and all elements of S are sorted.
  
  dd = packing_set_test(E, S) 
  first, second, third = S[1], S[2], S[3]
  required_states = [
                     [[first, second, third], [third]],
                     [[first, second, third], [second], [third]],
                     [[second, third], [second], [third]],
                     [[second, third], [second], [third], []],
                     [[third], []],
                     [[]],
                    ]
  for idx in 1:length(dd.variables)
    states = [state["uncovered_sets"] for state in dd.states[idx]]
    @test length(required_states[idx]) == length(states) && all([i in states for i in required_states[idx]])
  end
end

function indepedent_set_test(V, E)
  #=
  Building DD
  =#
  # parameters
  empty_domain = -1                            # any value is valid, since the domain 0 always will be available (check get_candidates function definition), this happens because the problem modeling guarantees this.
  domain = Array{Any}([0, 1])                  # 0, if we do not select the node, and 1 otherwise
  variables = Array{Any}(V)                    # the nodes
  initial_state = State("available_nodes" => V) # all nodes

  # callbacks
  get_next_state(idx_variable::Int64, state::DecisionDiagram.State, candidate::Any, dd::DecisionDiagram.Instance) = State("available_nodes" => filter(i -> dd.variables[idx_variable] != i && (candidate == 0 || !in((dd.variables[idx_variable], i), E) && !in((i, dd.variables[idx_variable]), E)), state["available_nodes"]))

  get_candidates(idx_variable::Int64, state::DecisionDiagram.State, dd::DecisionDiagram.Instance) = dd.variables[idx_variable] in state["available_nodes"] ? [0, 1] : [0] 

  not_redundant_state(idx_variable::Int64, state::State, states::Array{State}) = true

  dd = DecisionDiagram.Instance(empty_domain, domain, variables, initial_state, get_candidates, get_next_state, not_redundant_state)
  dd.LOG = true
  DecisionDiagram.run(dd)

  return dd
end

function indepedent_set_test_case_1()
  V, E = 1:5, [(1, 2), (1, 3), (2, 3), (2, 4), (3, 4), (4, 5)] # undirected graph
  
  dd = indepedent_set_test(V, E) 
  required_states = [
                     [[2, 3, 4, 5], [4, 5]],
                     [[5], [3, 4, 5], [4, 5]],
                     [[5], [4, 5]],
                     [[5], []],
                     [[]],
                    ]
  for idx in 1:length(dd.variables)
    states = [state["available_nodes"] for state in dd.states[idx]]
    @test length(required_states[idx]) == length(states) && all([i in states for i in required_states[idx]])
  end
end

function makespan_test(times)
  positions = jobs = 1:length(times)
  #=
  Building DD
  =#
  # parameters
  empty_domain = -1                             # any value is valid, since the domain 0 always will be available (check get_candidates function definition), this happens because the problem modeling guarantees this.
  domain = Array{Any}(jobs)                     # the job to be placed in the in the i-th position, for all i in variables
  variables = Array{Any}(positions)             # the available positions
  initial_state = State("allocated_jobs" => []) # initially we have an 

  # callbacks
  get_next_state(idx_variable::Int64, state::DecisionDiagram.State, candidate::Any, dd::DecisionDiagram.Instance) = State("allocated_jobs" => union(state["allocated_jobs"], [candidate])) 

  get_candidates(idx_variable::Int64, state::DecisionDiagram.State, dd::DecisionDiagram.Instance) = filter(job -> !in(job, state["allocated_jobs"]), dd.domain)

  not_redundant_state(idx_variable::Int64, state::State, states::Array{State}) = all([Set{Int64}(state["allocated_jobs"]) != Set{Int64}(state′["allocated_jobs"]) for state′ in states])

  dd = DecisionDiagram.Instance(empty_domain, domain, variables, initial_state, get_candidates, get_next_state, not_redundant_state)
  #dd.LOG = true
  DecisionDiagram.run(dd)

  return dd
end

function makespan_test_case_1()
  times = [[4, 5, 9], [3, 7, 8], [1, 2, 10]]
  
  dd = makespan_test(times) 
  required_states = [
                     [[1], [2], [3]],
                     [[1, 2], [1, 3], [2, 3]],
                     [[1, 2, 3]],
                    ]
  for idx in 1:length(dd.variables)
    states = [state["allocated_jobs"] for state in dd.states[idx]]
    @test length(required_states[idx]) == length(states) && all([i in states for i in required_states[idx]])
  end
end

function test()
  knapsack_test_case_1()
  knapsack_test_case_2()
  set_covering_test_case_1()
  packing_set_test_case_1()
  indepedent_set_test_case_1()
  makespan_test_case_1()
end

end
