module DecisionDiagram

export State

const State = Dict{Any, Any}

mutable struct Instance 
  empty_domain::Any
  domain::Vector{Any}
  variables::Vector{Any}
  initial_state::State
  width_limit::Int
  get_candidates
  get_next_state
  not_redundant_state
  compact_width 
  states::Vector{Vector{State}} # idx - variable index
  LOG::Bool
  Instance(empty_domain::Any, domain::Vector{Any}, variables::Vector{Any}, initial_state::State, get_candidates, get_next_state, not_redundant_state) = new(empty_domain, domain, variables, initial_state, typemax(Int), get_candidates, get_next_state, not_redundant_state, nothing, Vector{Vector{State}}([Vector{State, 1}() for variable in variables]), false)
end

function forward(idx_variable::Int, state::State, dd::Instance)
  candidates = dd.get_candidates(idx_variable, state, dd)  
  isempty(candidates) && push!(candidates, dd.empty_domain)
  dd.LOG && println("  - In $state we have the candidates $candidates")
  for candidate in candidates
    next_state = dd.get_next_state(idx_variable, state, candidate, dd)
    add_state(idx_variable, next_state, dd) 
    dd.LOG && println("   > Candidate $candidate generated state $next_state")
  end
end

function add_state(idx_variable::Int, state::State, dd::Instance)
  states = dd.states[idx_variable]
  (!in(state, states) && dd.not_redundant_state(idx_variable, state, states)) && push!(states, state)
end

function compact_width(idx_variable::Int, dd::Instance)
  n = length(dd.states[idx_variable])
  if n > dd.width_limit
    dd.compact_width == nothing && error("The index $idx_variable (with $n states) reached the width limit of $(dd.width_limit) and the `compact_width` function was not defined in the DecisionDiagram")
    dd.LOG && println(" * COMPACTING WIDTH FROM $n to $(dd.width_limit).")
    dd.compact_width(idx_variable)
    dd.LOG && println(" * NEW STATES ARE THESE $(dd.states[idx_variable]).")
  end
end

function run(dd::Instance)
  dd.LOG && println(" * Forwarding at the root.")
  forward(1, dd.initial_state, dd) # initialise in the first variable
  dd.LOG && println(" * First depth has the states $(dd.states[1]).")
  compact_width(1, dd)
  n = length(dd.variables)
  for idx in 1:n - 1 
    dd.LOG && println(" * Forwarding to index $(idx + 1) (variable $(dd.variables[idx + 1])):")
    states = dd.states[idx]
    isempty(states) && error("The instance is infeasible")
    [forward(idx + 1, state, dd) for state in states]
    dd.LOG && println(" * Depth $(idx + 1) has the states $(dd.states[idx + 1]).")
    compact_width(idx + 1, dd)
  end
  isempty(dd.states[n]) && error("The instance is infeasible")
end

end

