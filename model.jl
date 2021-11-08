using CPLEX
using JuMP
using BenchmarkTools
include("SparseMaxFlowMinCut.jl")

function bfs(x_val, i)
  δ⁺′(i), S, q = [a for (a, v) in x_val if a[1] == i && v > 0.001], Set{Int64}([i]), [i]
  while !isempty(q)
    curr = popfirst!(q)
    A′ = [(i, j) for (i, j) in δ⁺′(curr) if j in S]
    [push!(S, next) for (curr, next) in A′]
    [push!(q, next) for (curr, next) in A′ if !in(next, S)]
  end
  return S
end

function bfs_sets(x_val, nodes, source::Int64)
  sets = []
  for i in nodes
    S = bfs(x_val, i, source)
    !in(depot, S) && length(S) > 1 && push!(sets, S)
  end
  return sets
end

function gh_sets(x_val, source::Int64, targets, n::Int64)
  g, M, sets = SparseMaxFlowMinCut.ArcFlow[], 100000, Set{Set{Int64}}()
  # mouting graph
  [push!(g, SparseMaxFlowMinCut.ArcFlow(i, j, trunc(floor(value, digits=5) * M))) for (a, value) in x_val if value > 0.0001]
  # get subsets
  for i in nodes 
    maxFlow, flows, set = SparseMaxFlowMinCut.find_maxflow_mincut(SparseMaxFlowMinCut.Graph(n, g), source, i)
    (maxFlow / M) < (2 - 0.001) && push!(sets, Set{Int64}([i for i in 1:n if cut[i] == 1]))
  end
  return sets
end

function build_model_sbrp(data::SBRPData, app::Dict{String,Any})
  B, A, T, V, depot, Vb, gh_cuts, bfs_cuts = data.B, data.D.A, data.T, data.D.V, data.depot, [i for i in b for b in data.B], [], []
  # Formulation
  model = Model(CPLEX.Optimizer())
  MOI.set(model, MOI.NumberOfThreads(), 1)
  @variable(model, x[a in A], Int)
  @objective(model, Min, sum(time(data, a) * x[a] for a in A))
  @constraint(model, degree[i in V], sum(x[a] for a in δ⁻(A, i)) == sum(x[a] for a in δ⁺(A, i)))
  @constraint(sbrp.formulation, block[block in B], sum(sum(x[a] for a in δ⁺(A, i)) for i in block) >= 1)
  @constraint(sbrp.formulation, sum(x[a] for a in A) <= T - sum(time_block(data, block) for block in B))
  # connectivity
  # lazy
  function bfs_callback(cb_data)
    callback_node_status(cb_data, model) != MOI.CALLBACK_NODE_STATUS_INTEGER && return 
    x_val = callback_value(cb_data, x)
    # bfs
    sets = bfs_sets(x_val, Vb, depot)
    for S in sets
      push!(bfs_cuts, S)
      # add ineq
      [MOI.submit(model, MOI.LazyConstraint(cb_data), @build_constraint(sum(x[a] for a in δ⁺(A, S)) >= x[(i, j)])) for (i, j) in A if i in S && j in S]
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
      [MOI.submit(model, MOI.UserCut(cb_data), @build_constraint(sum(x[a′] for a′ in δ⁺(A, S)) >= x[a])) for a in δ⁺(S)]
    end
  end
  MOI.set(model, CPLEX.CallbackFunction(), gh_callback)
  return (model, x)
end

function build_model_max_profit_sbrp(data::SBRPData, app::Dict{String,Any})
  B, A, T, V, depot, Vb, gh_cuts, bfs_cuts = data.B, Set{Int64}(data.D.A), data.T, data.D.V, data.depot, Set{Int64}([i for i in b for b in data.B]), [], []
  δ⁺(S) = [(i, j) for (i, j) in A if i in S && !in(j, S)]
  get_block(i) = find(block -> i in block, B)
  # Formulation
  model = Model(CPLEX.Optimizer)
  @variable(model, x[a in A], Int)
  @variable(model, y[block in B], Int)
  @objective(model, Min, sum(time(data, a) * x[a] for a in A))
  @constraint(model, degree[i in V], sum(x[a] for a in δ⁻(A, i)) == sum(x[a] for a in δ⁺(A, i)))
  @constraint(sbrp.formulation, block[block in B], sum(sum(x[a] for a in δ⁺(A, i)) for i in block) >= y[b])
  @constraint(sbrp.formulation, sum(x[a] for a in A) + sum(y[block] * time_block(data, block) for block in B) <= T)
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
