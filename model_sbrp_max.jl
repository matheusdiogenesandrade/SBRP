module ModelSBRPMax

using ..Model
using ...Data
using ...Data.SBRP
using CPLEX
using JuMP

export build_model_sbrp_max

blocks(B, S) = Array{Array{Int64, 1}, 1}([b for b in B if all([i in S for i in b])])

function add_subtour_ineqs(model, x, y, sets, A::Array{Tuple{Int64, Int64}})
  [@constraint(model, sum(x[a′] for a′ in δ⁺(A, S)) >= y[block]) for (S, block) in sets]
end

function get_max_flow_min_cut_cuts(data::SBRPData, model, x, y)
  A, B, depot, sets = data.D.A, data.B, data.depot, Set{Tuple{Set{Int64}, Array{Int64, 1}}}()
  n = max([max(a...) for a in A]...)
  while true
    optimize!(model)
    !in(termination_status(model), [MOI.OPTIMAL, MOI.TIME_LIMIT]) && error("The model could not be solved")
    # create dummy nodes
    Nb = Dict{Array{Int64, 1}, Int64}(B[i] => n + i for i in 1:length(B))
    ids = Set{Int64}(collect(values(Nb)))
    n′ = n + length(B)
    # get values
    y_val = Dict{Array{Int64, 1}, Float64}(block => value(y[block]) for block in B)
    # get subsets
    sets′ = Set{Set{Int64}}()
    g, M, sets″ = SparseMaxFlowMinCut.ArcFlow[], 100000, Set{Set{Int64}}()
    # mouting graph
    [push!(g, SparseMaxFlowMinCut.ArcFlow(a[1], a[2], trunc(floor(value(x[a]), digits=5) * M))) for a in A if value(x[a]) > 1e-3]
    [push!(g, SparseMaxFlowMinCut.ArcFlow(i, j, 1 * M)) for (b, i) in Nb for j in b]
    # get subsets
    for source in ids 
      maxFlow, flows, set = SparseMaxFlowMinCut.find_maxflow_mincut(SparseMaxFlowMinCut.Graph(n′, g), source, depot)
      flow = maxFlow / M
      # base case
      (flow >= 1 - 1e-3 && set[depot] == 0) && continue
      # get set
      S = Set{Int64}([i for i in 1:n′ if set[i] == 1])
      # get blocks inside
      blocks = blocks(S)
      # check validity
      [push!(sets′, (setdiff(S, ids), block)) for block in blocks if y_val[block] > flow]
    end
    # base case
    isempty(sets′) && break
#    println(sets′)
    # store sets″
    union!(sets, sets′)
    # add ineqs
    add_subtour_ineqs(model, x, y, sets′, A)
  end
  return sets
end

function build_model_sbrp_max(data::SBRPData, app::Dict{String,Any})
  B, A, T, V, depot, profits, Vb, info = data.B, data.D.A, data.T, 1:length(data.D.V), data.depot, data.profits, Set{Int64}([i for b in data.B for i in b]), Dict{String, Any}("lazyCuts" => 0)
  function add_basic_constraints(model)
    @objective(model, Max, sum(data.profits[b] * y[b] for b in B))
    @constraint(model, degree[i in V], sum(x[a] for a in δ⁻(A, i)) == sum(x[a] for a in δ⁺(A, i)))
    @constraint(model, sum(x[a] for a in δ⁺(A, depot)) <= 1)
    @constraint(model, block[block in B, i in block], sum(x[a] for a in δ⁺(A, i)) >= y[block])
    @constraint(model, block1[block in B], sum(sum(x[a] for a in δ⁺(A, i)) for i in block) >= y[block])
    @constraint(model, block2[block in B], sum(x[a] for a in δ⁺(A, block)) >= y[block])
#    @constraint(model, sum(x[a] for a in A) <= T - sum(time_block(data, block) for block in B))
  end
  # Formulation
  # frac model - get max-flow/min-cut cuts
  model = direct_model(CPLEX.Optimizer())
  set_silent(model)
  @variable(model, x[a in A], lower_bound = 0, upper_bound = 1)
#  @variable(model, y[b in B], lower_bound = 0, upper_bound=1)
  @variable(model, y[b in B], Bin)
  add_basic_constraints(model)
  info["maxFlowCutsTime"] = @elapsed sets = get_max_flow_min_cut_cuts(data, model, x, y)
  info["maxFlowCuts"], info["rootLP"] = length(sets), objective_value(model)
  # integer model
  model = direct_model(CPLEX.Optimizer())
  set_silent(model)
  MOI.set(model, MOI.NumberOfThreads(), 1)
  @variable(model, x[a in A], Bin)
#  @variable(model, y[b in B], lower_bound = 0, upper_bound=1)
  @variable(model, y[b in B], lower_bound = 0, upper_bound = 1)
  add_basic_constraints(model)
  add_subtour_ineqs(model, x, y, sets, A)
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
    y_val = callback_value.(Ref(cb_data), y)
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
      # add ineq
      blocks_ = blocks(B, S′)
      for block in blocks_
        abs(y_val[block] - length(δ⁺(A′, S′))) < 1e-3 && continue
        println(collect(S′), " ", block)
        MOI.submit(model, MOI.LazyConstraint(cb_data), @build_constraint(sum(x[a] for a in δ⁺(A, S′)) >= y[block]))
        info["lazyCuts"] = info["lazyCuts"] + 1
      end
    end
  end
  MOI.set(model, CPLEX.CallbackFunction(), bfs_callback)
  return (model, x, y, info)
end

end
