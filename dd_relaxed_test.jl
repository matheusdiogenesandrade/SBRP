module RelaxedDecisionDiagramTests

include("dd.jl")

using .DecisionDiagram
using Test

function indepedent_set_test(V, E)
  #=
  Building DD
  =#
  # parameters
  empty_domain = -1                            # any value is valid, since the domain 0 always will be available (check get_candidates function definition), this happens because the problem modeling guarantees this.
  domain = Vector{Any}([0, 1])                  # 0, if we do not select the node, and 1 otherwise
  variables = Vector{Any}(V)                    # the nodes
  initial_state = State("available_nodes" => V) # all nodes

  # callbacks
  get_next_state(idx_variable::Int, state::DecisionDiagram.State, candidate::Any, dd::DecisionDiagram.Instance) = State("available_nodes" => filter(i -> dd.variables[idx_variable] != i && (candidate == 0 || !in((dd.variables[idx_variable], i), E) && !in((i, dd.variables[idx_variable]), E)), state["available_nodes"]))

  get_candidates(idx_variable::Int, state::DecisionDiagram.State, dd::DecisionDiagram.Instance) = dd.variables[idx_variable] in state["available_nodes"] ? [0, 1] : [0] 

  not_redundant_state(idx_variable::Int, state::State, states::Vector{State}) = true

  dd = DecisionDiagram.Instance(empty_domain, domain, variables, initial_state, get_candidates, get_next_state, not_redundant_state)
  dd.LOG = true

  # settings for relaxed DD
  dd.width_limit = 2 
  function compact_width(idx_variable::Int) 
    states = dd.states[idx_variable]
    states_to_compact = states[1:length(states) - dd.width_limit + 1]
    setdiff!(states, states_to_compact)
    compacted_state = State("available_nodes" => union(i for state in states_to_compact for i in state["available_nodes"]))
    push!(states, compacted_state)
  end
  dd.compact_width = compact_width

  # run
  
  DecisionDiagram.run(dd)

  return dd
end

function indepedent_set_test_case_1()
  V, E = 1:5, [(1, 2), (1, 3), (2, 3), (2, 4), (3, 4), (4, 5)] # undirected graph
  
  dd = indepedent_set_test(V, E) 
  required_states = [
                     [[2, 3, 4, 5], [4, 5]],
                     [[3, 4, 5], [4, 5]],
                     [[5], [4, 5]],
                     [[5], []],
                     [[]],
                    ]
  for idx in 1:length(dd.variables)
    states = [state["available_nodes"] for state in dd.states[idx]]
    @test length(required_states[idx]) == length(states) && all([i in states for i in required_states[idx]])
  end
end

function test()
  indepedent_set_test_case_1()
end

end
