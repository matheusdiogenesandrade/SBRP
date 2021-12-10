module Model

using ..Data
using ..Data.SBRP
using CPLEX
using JuMP
using BenchmarkTools

export ModelSBRP, add_subtour_ineqs, SparseMaxFlowMinCut, EPS, lazy_separation, get_max_flow_min_cut_cuts, bfs_sets, blocks

include("SparseMaxFlowMinCut.jl")

EPS = 1e-3
blocks(B, S) = Array{Array{Int64, 1}, 1}([b for b in B if any([i in S for i in b])])
add_subtour_ineqs(model, x, y, sets::Set{Tuple{Set{Int64}, Array{Int64, 1}}}, A::Array{Tuple{Int64, Int64}, 1}) = [@constraint(model, sum(x[a] for a in δ⁺(A, S)) >= y[block]) for (S, block) in sets]

function bfs(A, i)
  S, q = Set{Int64}([i]), [i]
  δ⁺′(j) = [a for a in A if a[1] == j && !in(a[2], S)]
  while !isempty(q)
    curr = popfirst!(q)
    A′ = δ⁺′(curr)
    [(push!(S, next), push!(q, next)) for (curr, next) in A′]
  end
  return S
end

function bfs_sets(A, nodes, source::Int64)
  sets = Set{Set{Int64}}()
  for i in nodes
    S = bfs(A, i)
    (!in(source, S) && length(S) > 1) && push!(sets, S)
  end
  return sets
end

function lazy_separation(data::SBRPData, Vb, info, model, x, y, cb_data::CPLEX.CallbackContext, context_id::Clong)
  # preliminary checkings
  context_id != CPX_CALLBACKCONTEXT_CANDIDATE && return
  ispoint_p = Ref{Cint}()
  (CPXcallbackcandidateispoint(cb_data, ispoint_p) != 0 || ispoint_p[] == 0) && return 
  # get values
  CPLEX.load_callback_variable_primal(cb_data, context_id)
  x_val, y_val = callback_value.(Ref(cb_data), x), callback_value.(Ref(cb_data), y)
  # bfs
  A′ = [a for a in data.D.A if x_val[a] > 0.5]
  sets = bfs_sets(A′, Vb, data.depot)
  # get valid components
  sets′ = Set{Tuple{Set{Int64}, Array{Int64, 1}}}((S, block) for S in sets for block in blocks(data.B, S) if y_val[block] > length(δ⁺(A′, S)) + 0.5)
  # add ineqs
  [MOI.submit(model, MOI.LazyConstraint(cb_data), @build_constraint(sum(x[a] for a in δ⁺(data.D.A, S)) >= y[block])) for (S, block) in sets′]
  # update info
  info["lazyCuts"] = info["lazyCuts"] + length(sets′)
  #
#  [println(collect(S), block) for (S, block) in sets′]
end

function get_max_flow_min_cut_cuts(data::SBRPData, model, x, y, info::Dict{String, Any})
  A, B, depot, sets, Vb = data.D.A, data.B, data.depot, Set{Tuple{Set{Int64}, Array{Int64, 1}}}(), Set{Int64}([i for b in data.B for i in b])
  n, iteration = max(Set{Int64}(vcat([i for (i, j) in A], [j for (i, j) in A]))...), 1
  while true
    time = @elapsed optimize!(model)
    !in(termination_status(model), [MOI.OPTIMAL, MOI.TIME_LIMIT]) && error("The model could not be solved")
    # get values
    y_val = Dict{Array{Int64, 1}, Float64}(block => value(y[block]) for block in B)
    # get subsets
    g, M, sets′ = SparseMaxFlowMinCut.ArcFlow[], 100000, Set{Tuple{Set{Int64}, Array{Int64, 1}}}()
    # mouting graph
    [push!(g, SparseMaxFlowMinCut.ArcFlow(a[1], a[2], trunc(floor(value(x[a]), digits=5) * M))) for a in A if value(x[a]) > EPS]
    # get subsets
    for source in Vb 
      maxFlow, flows, set = SparseMaxFlowMinCut.find_maxflow_mincut(SparseMaxFlowMinCut.Graph(n, g), source, depot)
      flow = maxFlow / M
      # base case
      (flow >= 1 - 1e-2 && set[depot] == 0) && continue
      # get set
      S = Set{Int64}([i for i in 1:n if set[i] == 1])
      # get blocks
      [push!(sets′, (S, block)) for block in blocks(B, S) if y_val[block] > flow + 1e-2]
    end
    # base case
    isempty(sets′) && break
#    [println(S, block) for (S, block) in sets′]
    # update info
    info["iteration_" * string(iteration) * "_time"] = time
    info["iteration_" * string(iteration) * "_cuts"] = length(sets′)
    # increment iteration
    iteration += 1
    # store sets″
    union!(sets, sets′)
    # add ineqs
    add_subtour_ineqs(model, x, y, sets′, A)
  end
  return sets
end

include("model_sbrp_max.jl")

include("model_sbrp_max_complete.jl")

end
