module Model

using ..Data
using ..Data.SBRP
using CPLEX
using JuMP
using BenchmarkTools

export ModelSBRP, get_max_flow_min_cut_cuts, bfs_sets, add_subtour_ineqs, SparseMaxFlowMinCut

include("SparseMaxFlowMinCut.jl")
#include("data.jl")

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

function gh_sets(x_val, sink::Int64, sources, n::Int64)
  #=
  for (a, v) in x_val
    println(a, " ", v)
  end
  =#
  g, M, sets = SparseMaxFlowMinCut.ArcFlow[], 100000, Set{Set{Int64}}()
  # mouting graph
#  [push!(g, SparseMaxFlowMinCut.ArcFlow(a[1], a[2], trunc(floor(value, digits=5) * M))) for (a, value) in x_val if value > 1e-3]
  [push!(g, SparseMaxFlowMinCut.ArcFlow(a[1], a[2], trunc(floor(v, digits=5) * M))) for (a, v) in x_val]
  [push!(g, SparseMaxFlowMinCut.ArcFlow(a[2], a[1], trunc(floor(v, digits=5) * M))) for (a, v) in x_val if !in(a[1], sources)]
  # get subsets
  for source in sources 
    maxFlow, flows, set = SparseMaxFlowMinCut.find_maxflow_mincut(SparseMaxFlowMinCut.Graph(n, g), source, sink)
#    ((maxFlow / M) < (2 - 1e-3) && set[sink] == 0) && println(source, " to ", sink, " = ", maxFlow / M, " ", set[sink], " ", Set{Int64}([i for i in 1:n if set[i] == 1]))
    ((maxFlow / M) < (2 - 1e-3) && set[sink] == 0) && push!(sets, Set{Int64}([i for i in 1:n if set[i] == 1]))
  end
  return sets
end

function add_subtour_ineqs(model, x, sets, A::Array{Tuple{Int64, Int64}})
  [@constraint(model, sum(x[a′] for a′ in δ⁺(A, S)) >= 1) for S in sets]
end

function get_max_flow_min_cut_cuts(data::SBRPData, model, x)
  A, B, depot, sets = data.D.A, data.B, data.depot, Set{Set{Int64}}()
  n = max([max(a...) for a in A]...)
  while true
    optimize!(model)
    !in(termination_status(model), [MOI.OPTIMAL, MOI.TIME_LIMIT]) && error("The model could not be solved")
    # create dummy nodes
    Nb = Dict{Array{Int64, 1}, Int64}(B[i] => n + i for i in 1:length(B))
    ids = Set{Int64}(collect(values(Nb)))
    # get value
    x_val = merge(
                  Dict{Tuple{Int64, Int64}, Float64}(a      => value(x[a]) for a in A if value(x[a]) > 1e-3),
                  Dict{Tuple{Int64, Int64}, Float64}((i, j) => 2.0         for (b, i) in Nb for j in b)
                 )
    # get subsets
    sets′, sets″ = gh_sets(x_val, depot, ids, n + length(B)), []
    for S in sets′
      S′ = setdiff(S, ids)
      !(isempty(S′) || in(S′, sets)) && push!(sets″, S′)
    end
    # base case
    isempty(sets″) && break
#    println(sets″)
    # store sets″
    [push!(sets, S) for S in sets″]
    # add ineqs
    add_subtour_ineqs(model, x, sets″, A)
  end
  return sets
end

#=
function build_model_atsp(V::Array{Int64}, A::Array{Tuple{Int64, Int64}}, costs::Dict{Tuple{Int64, Int64}, Float64}, depot::Int64)
  Vₙ = [i for i in V if i != depot]
  model = direct_model(CPLEX.Optimizer())
  set_silent(model)
  @variable(model, x[a in A], Bin)
  @variable(model, y[a in A], lower_bound=0, upper_bound=length(V))
  @objective(model, Min, sum(costs[a] * x[a] for a in A))
  @constraint(model, degree[i in V], sum(x[a] for a in δ⁻(A, i)) == sum(x[a] for a in δ⁺(A, i)))
  @constraint(model, degree_once[i in V], sum(x[a] for a in δ⁻(A, i)) == 1)
  @constraint(model, mtz[i in Vₙ], sum(y[a] for a in δ⁺(A, i)) == sum(y[a] for a in δ⁻(A, i)) + sum(x[a] for a in δ⁺(A, i)))
  @constraint(model, mtz1[a in A], y[a] >= x[a])
  @constraint(model, mtz2[a in A], x[a] * length(V) >= y[a])
  return (model, x, y)
end
=#

#=
function build_model_max_profit_sbrp(data::SBRPData, app::Dict{String,Any})
B, A, T, V, depot, Vb, gh_cuts, bfs_cuts = data.B, Set{Int64}(data.D.A), data.T, data.D.V, data.depot, Set{Int64}([i for i in b for b in data.B]), [], []
δ⁺(S) = [(i, j) for (i, j) in A if i in S && !in(j, S)]
get_block(i) = find(block -> i in block, B)
# Formulation
model = Model(CPLEX.Optimizer)
@variable(model, x[a in A], Int, lower_bound = 0)
@variable(model, y[block in B], Bin)
@objective(model, Min, sum(time(data, a) * x[a] for a in A))
@constraint(model, degree[i in V], sum(x[a] for a in δ⁻(A, i)) == sum(x[a] for a in δ⁺(A, i)))
@constraint(sbrp.formulation, block[block in B], sum(sum(x[a] for a in δ⁺(A, i)) for i in block) >= y[b])
@constraint(sbrp.formulation, sum(x[a] for a in A) + sum(y[block] * time_block(data, block) for block in B) <= T)
@constraint(sbrp.formulation, sum(x[a] for a in δ⁺(A, depot)) <= 1)
@constraint(sbrp.formulation, [(block, i, a) in [(block, i, a) for block in B for i in b for a in δ⁺(A, i)]], y[block] >= x[a])
# connectivity
# lazy
function bfs_callback(cb_data)
callback_node_status(cb_data, model) != MOI.CALLBACK_NODE_STATUS_INTEGER && continue
x_val = callback_value(cb_data, x), 
for S in bfs_sets(x_val, Vb, depot)
push!(bfs_cuts, S)
[MOI.submit(model, MOI.LazyConstraint(cb_data), @build_constraint(sum(x[a] for a in δ⁺(A, S)) >= term)) for term in [x[get_block(i)] for i in intersect(S, Vb)]]
end
end
MOI.set(model, MOI.LazyConstraintCallback(), bfs_callback)
# user cut
function gh_callback(cb_data::CPLEX.CallbackContext, context_id::Clong)
# if it is in the root
depth = Ref{Cint}()
#ret = 
CPXcallbackgetinfodbl(cb_data, CPX_CALLBACK_INFO_NODE_DEPTH, depth)
depth > 0 && return
# get valus
CPLEX.load_callback_variable_primal(cb_data, context_id)
x_val = callback_value(cb_data, x)
# gh
sets = gh_sets(x_val, depot, Vb, length(V))
for S in sets
S in gh_cuts && continue
push!(gh_cuts, S)
[MOI.submit(model, MOI.UserCut(cb_data), @build_constraint(sum(x[a] for a in δ⁺(A, S)) >= y[block])) for block in Set{Array{Int64}}([get_block(i) for i in intersect(S, Vb)])]
end
end
MOI.set(model, CPLEX.CallbackFunction(), gh_callback)
return (model, x, y)
end
=#

include("model_sbrp.jl")

include("model_sbrp_complete.jl")

include("model_tsp.jl")

include("model_sbrp_max.jl")

end
