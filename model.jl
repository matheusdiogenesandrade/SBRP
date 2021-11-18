using CPLEX
using JuMP
using BenchmarkTools
include("SparseMaxFlowMinCut.jl")
include("data.jl")

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
  [push!(g, SparseMaxFlowMinCut.ArcFlow(a[2], a[1], trunc(floor(v, digits=5) * M))) for (a, v) in x_val]
  # get subsets
  for source in sources 
    maxFlow, flows, set = SparseMaxFlowMinCut.find_maxflow_mincut(SparseMaxFlowMinCut.Graph(n, g), source, sink)
#    println(source, " to ", sink, " = ", maxFlow / M, " ", set[sink])
    ((maxFlow / M) < (2 - 1e-3) && set[sink] == 0) && push!(sets, Set{Int64}([i for i in 1:n if set[i] == 1]))
  end
  return sets
end

function get_max_flow_min_cut_cuts(data::SBRPData, model, x)
  A, B, depot, sets, n = data.D.A, data.B, data.depot, Set{Set{Int64}}(), length(data.D.V)
  while true
    !solve(model) && error("The model could not be solved")
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
    println(sets″)
    # store sets″
    [push!(sets, S) for S in sets″]
    # add ineqs
    add_subtour_ineqs(model, x, sets″, A)
  end
  return sets
end

function add_subtour_ineqs(model, x, sets, A::Array{Tuple{Int64, Int64}})
  [@constraint(model, sum(x[a′] for a′ in δ⁺(A, S)) >= 1) for S in sets]
end

function build_model_sbrp(data::SBRPData, app::Dict{String,Any})
  B, A, T, V, depot, Vb, bfs_cuts = data.B, data.D.A, data.T, 1:length(data.D.V), data.depot, Set{Int64}([i for b in data.B for i in b]), [], Set{Set{Int64}}()
  function add_basic_constraints(model)
    @objective(model, Min, sum(time(data, a) * x[a] for a in A))
    @constraint(model, degree[i in V], sum(x[a] for a in δ⁻(A, i)) == sum(x[a] for a in δ⁺(A, i)))
    @constraint(model, block[block in B], sum(sum(x[a] for a in δ⁺(A, i)) for i in block) >= 1)
    @constraint(model, block1[block in B], sum(x[a] for a in δ⁺(A, block))  >= 1)
    @constraint(model, sum(x[a] for a in δ⁺(A, depot)) <= 1)
#    @constraint(model, sum(x[a] for a in A) <= T - sum(time_block(data, block) for block in B))
  end
  # Formulation
  # frac model - get max-flow/min-cut cuts
  model = direct_model(CPLEX.Optimizer())
  set_silent(model)
  @variable(model, x[a in A], lower_bound = 0)
  add_basic_constraints(model)
  sets = get_max_flow_min_cut_cuts(data, model, x)
  println("# Max-flow-min-cuts: $(length(sets))")
  # integer model
  model = direct_model(CPLEX.Optimizer())
  set_silent(model)
  MOI.set(model, MOI.NumberOfThreads(), 1)
  @variable(model, x[a in A], Int, lower_bound = 0)
  add_basic_constraints(model)
  add_subtour_ineqs(model, x, sets, A)
  # connectivity
  # lazy
  function bfs_callback(cb_data::CPLEX.CallbackContext, context_id::Clong)
    # preliminary checkings
    context_id != CPX_CALLBACKCONTEXT_CANDIDATE && return
    ispoint_p = Ref{Cint}()
    (CPXcallbackcandidateispoint(cb_data, ispoint_p) != 0 || ispoint_p[] == 0) && return 
    # get values
    CPLEX.load_callback_variable_primal(cb_data, context_id)
    x_val = callback_value.(Ref(cb_data), x)
    # create dummy nodes
    n = length(V)
    Nb = Dict{Array{Int64, 1}, Int64}(B[i] => n + i for i in 1:length(B))
    ids = Set{Int64}(collect(values(Nb)))
    # bfs
    A′ = vcat([a for a in A if x_val[a] > 0.5], [(i, j) for (b, i) in Nb for j in b])
    sets = bfs_sets(A′, ids, depot)
#    !isempty(sets) && println("=========Lazy callback==========")
    for S in sets
      # remove dummy node
      S′ = setdiff(S, ids)
      # edge case
      isempty(S′) && continue
      # store
      push!(bfs_cuts, S′)
      println(collect(S′))
      # add ineq
      MOI.submit(model, MOI.LazyConstraint(cb_data), @build_constraint(sum(x[a] for a in δ⁺(A, S′)) >= 1))
    end
  end
  MOI.set(model, CPLEX.CallbackFunction(), bfs_callback)
  return (model, x)
end

function build_model_sbrp_complete(data::SBRPData, app::Dict{String,Any})
  B, A, T, depot, Vb, gh_cuts, bfs_cuts = data.B, data.D.A, data.T, data.depot, Set{Int64}([i for b in data.B for i in b]), [], []
  V = Set{Int64}(vcat([i for (i, j) in A], [j for (i, j) in A]))
  # Formulation
  model = direct_model(CPLEX.Optimizer())
  set_silent(model)
#  MOI.set(model, MOI.NumberOfThreads(), 1)
  @variable(model, x[a in A], Bin)
  @variable(model, y[a in A], lower_bound=0, upper_bound=length(V))
  @objective(model, Min, sum(time(data, a) * x[a] for a in A))
  @constraint(model, degree[i in V], sum(x[a] for a in δ⁻(A, i)) == sum(x[a] for a in δ⁺(A, i)))
  @constraint(model, sum(x[a] for a in δ⁺(A, depot)) <= 1)
  @constraint(model, mtz[i in Vb], sum(y[a] for a in δ⁺(A, i)) == sum(y[a] for a in δ⁻(A, i)) + sum(x[a] for a in δ⁺(A, i)))
  @constraint(model, mtz1[a in A], y[a] >= x[a])
  @constraint(model, mtz2[a in A], x[a] * length(V) >= y[a])
  @constraint(model, block[block in B], sum(sum(x[a] for a in δ⁺(A, i)) for i in block) == 1)
#  @constraint(model, sum(x[a] * time(data, a) for a in A) <= T - sum(time_block(data, block) for block in B))
  return (model, x, y)
end

function build_atsp_instance(data::SBRPData)
  # polynomial reduction 
  # SBRP attrs
  B, A, T, depot = data.B, data.D.A, data.T, data.depot
  V = Set{Int64}(vcat([i for (i, j) in A], [j for (i, j) in A]))
  # TSP attrs
  Vb, Vb′, n = Dict{Tuple{Int64, Array{Int64, 1}}, Int64}(), Dict{Tuple{Int64, Array{Int64, 1}}, Int64}(), max(collect(V)...) + 1
  # dummy nodes
  [(Vb[(i, b)] = n; n = n + 1) for b in B for i in b]
  [(Vb′[(i, b)] = n; n = n + 1) for b in B for i in b]
  # weights
  costs = merge!(
    Dict{Tuple{Int64, Int64}, Float64}((Vb[(i, b)],       Vb′[(i, b)])       => 0                       for b in B for i in b), # cycle arcs
    Dict{Tuple{Int64, Int64}, Float64}((Vb′[(b[i], b)],   Vb[(b[i + 1], b)]) => 0                       for b in B for i in 1:(length(b) - 1)), # cycle arcs
    Dict{Tuple{Int64, Int64}, Float64}((Vb′[(b[end], b)], Vb[(b[begin], b)]) => 0                       for b in B), # cycle arcs
    Dict{Tuple{Int64, Int64}, Float64}((Vb[(i, b)],       depot)             => 0                       for b in B for i in b), # depor arcs 
    Dict{Tuple{Int64, Int64}, Float64}((depot,            Vb′[(i, b)])       => 0                       for b in B for i in b), # depor arcs 
    Dict{Tuple{Int64, Int64}, Float64}((Vb[(i, b)],       Vb′[(j, b′)])      => time(data, (i, j)) for b in B for b′ in B for i in b for j in b′ if b != b′ && (i, j) in keys(data.D.distance)) # block arcs
  )
  # mip atsp
  Aᵗ = collect(keys(costs))
  Vᵗ = collect(Set{Int64}(vcat([i for (i, j) in Aᵗ], [j for (i, j) in Aᵗ])))
  return Vᵗ, Aᵗ, costs, Vb, Vb′
end

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
