module Model

using ..Symbols
using ..Data
using ..Data.SBRP
using ..NearestNeighborhoodHeuristic
using CPLEX
using JuMP
using BenchmarkTools

export ModelSBRP, add_subtour_ineqs, SparseMaxFlowMinCut, EPS, lazy_separation, get_max_flow_min_cut_cuts, bfs_sets, blocks, get_cc_cuts, add_cc_ineqs, add_intersection_cuts

include("SparseMaxFlowMinCut.jl")

EPS = 1e-4
blocks(B::Array{Array{Int64, 1}, 1}, S::Set{Int64}) = Array{Array{Int64, 1}, 1}([b for b in B if any([i in S for i in b])])
#blocks(B::Array{Array{Int64}}, i::Int64) = Array{Array{Int64, 1}, 1}([b for b in B if i in b])
add_cc_ineqs(model, x, A::Array{Tuple{Int64, Int64}, 1}, sets::Set{Int64}...) = [@constraint(model, sum(x[a] for a in δ⁺(A, S)) >= 2.0) for S in sets]
add_subtour_ineqs(model, x, y, sets::Set{Tuple{Set{Int64}, Int64, Array{Int64, 1}}}, A::Array{Tuple{Int64, Int64}, 1}) = [@constraint(model, sum(x[a] for a in δ⁺(A, S)) >= sum(x[(i, j)] for (i, j) in δ⁺(A, i) if j in S) + y[block] - 1) for (S, i, block) in sets]

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
  A = data.D.A
  A′ = [a for a in A if x_val[a] > 0.5]
  sets = bfs_sets(A′, Vb, data.depot)
  # get valid components
  sets′ = Set{Tuple{Set{Int64}, Int64, Array{Int64, 1}}}((S, i, block) for S in sets for block in blocks(data.B, S) for i in intersect(block, S) if y_val[block] + sum(x_val[(i, j)] for (i, j) in δ⁺(A, i) if j in S && x_val[(i, j)] > EPS) - 1 > length(δ⁺(A′, S)) + 0.5)
  # add ineqs
  [MOI.submit(model, MOI.LazyConstraint(cb_data), @build_constraint(sum(x[a] for a in δ⁺(A, S)) >= y[block] + sum(x[(i, j)] for (i, j) in δ⁺(A, i) if j in S) - 1)) for (S, i, block) in sets′]
  # update info
  info["lazyCuts"] = info["lazyCuts"] + length(sets′)
  #
  [println(collect(S), i, block) for (S, i, block) in sets′]
  flush(stdout)
end

function get_max_flow_min_cut_cuts(data::SBRPData, model, x, y, info::Dict{String, Any})
  A, B, depot, sets, Vb = data.D.A, data.B, data.depot, Set{Tuple{Set{Int64}, Int64, Array{Int64, 1}}}(), Set{Int64}([i for b in data.B for i in b])
  n, iteration = max(Set{Int64}(vcat([i for (i, j) in A], [j for (i, j) in A]))...), 1
  while true
    time = @elapsed optimize!(model)
    !in(termination_status(model), [MOI.OPTIMAL, MOI.TIME_LIMIT, MOI.ALMOST_INFEASIBLE]) && error("The model could not be solved")
    # get values
    y_val = Dict{Array{Int64, 1}, Float64}(block => value(y[block]) for block in B)
    x_val = Dict{Tuple{Int64, Int64}, Float64}(a => value(x[a]) for a in A)
    # get subsets
    g, M, sets′ = SparseMaxFlowMinCut.ArcFlow[], 100000, Set{Tuple{Set{Int64}, Int64, Array{Int64, 1}}}()
    # mouting graph
    [push!(g, SparseMaxFlowMinCut.ArcFlow(a[1], a[2], trunc(floor(x_val[a], digits=5) * M))) for a in A if x_val[a] > EPS]
    # get subsets
    for source in Vb 
      maxFlow, flows, set = SparseMaxFlowMinCut.find_maxflow_mincut(SparseMaxFlowMinCut.Graph(n, g), source, depot)
      flow = maxFlow / M
      # base case 1
      (flow >= 1 - 1e-2 || set[depot] == 1) && continue
      # get set
      S = Set{Int64}([i for i in 1:n if set[i] == 1])
      # base case 2
      length(S) <= 1 && continue
      # get blocks
      [push!(sets′, (S, i, block)) for block in blocks(B, S) for i in intersect(block, S) if y_val[block] + sum(x_val[(i, j)] for (i, j) in δ⁺(A, i) if j in S) - 1 > flow + 1e-2]
    end
    # base case
    isempty(sets′) && break
    println(length(sets′))
    flush(stdout)
#    [println(S, i, block) for (S, i, block) in sets′]
    # update infos
    info["iteration_" * string(iteration) * "_time"], info["iteration_" * string(iteration) * "_cuts"], iteration = time, length(sets′), iteration + 1
    # store sets″
    union!(sets, sets′)
    # add ineqs
    add_subtour_ineqs(model, x, y, sets′, A)
  end
  info["maxFlowCuts"], info["maxFlowLP"] = length(sets), objective_value(model)
  return sets
end

function get_cc_cuts(data::SBRPData, model, x, y, info::Dict{String, Any})
  # setup
  A, depot, B, sets, LB, iteration = data.D.A, data.depot, data.B, Set{Set{Int64}}(), NearestNeighborhoodHeuristic.solve(data)[1], 1
  V = Set{Int64}(vcat([i for (i, j) in A], [j for (i, j) in A]))
  while true
    time = @elapsed optimize!(model)
    !in(termination_status(model), [MOI.OPTIMAL, MOI.TIME_LIMIT]) && error("The model could not be solved")
    # get values
    x_val = Dict{Tuple{Int64, Int64}, Float64}(a => value(x[a]) for a in A)
    # model
    model_slave = direct_model(CPLEX.Optimizer())
    set_silent(model_slave)
    @variable(model_slave, z[a in A], Bin)
    @variable(model_slave, w[i in V], Bin)
    @variable(model_slave, y[block in B], Bin)
    @objective(model_slave, Min, sum(x_val[a] * z[a] for a in A))
    @constraint(model_slave, ub_z[(i, j) in A], w[i] + w[j] >= z[(i, j)])
    @constraint(model_slave, ub1_z[(i, j) in A], 2 - w[i] - w[j] >= z[(i, j)])
    @constraint(model_slave, lb_z[(i, j) in A], z[(i, j)] >= w[i] - w[j])
    @constraint(model_slave, lb1_z[(i, j) in A], z[(i, j)] >= w[j] - w[i])
    @constraint(model_slave, w[depot] == 1)
    @constraint(model_slave, sum(w[i] for i in V) >= 2)
    @constraint(model_slave, ub_y[block in B], sum(w[i] for i in block) >= y[block])
    @constraint(model_slave, sum(y[block] * data.profits[block] for block in B) <= LB - 0.5)
    @constraint(model_slave, sum(y[block] * time_block(data, block) for block in B) <= data.T)
    # solve
    optimize!(model_slave)
    !in(termination_status(model_slave), [MOI.OPTIMAL, MOI.TIME_LIMIT]) && error("The model could not be solved")
    # base case
    objective_value(model_slave) >= 2.0 - EPS && break
    # update infos
    info["cc_iteration_" * string(iteration) * "_time"], iteration = time, iteration + 1
    # store cut
    S = Set{Int64}(i for i in V if value(w[i]) > 0.5)
    push!(sets, S)
    add_cc_ineqs(model, x, A, [S]...)

    println(objective_value(model_slave), ": ", collect(S), ", ", length(S), ", ", δ⁺(A, S))
    flush(stdout)
  end
  info["cc_cuts"] = length(sets)
  return sets
end

function add_intersection_cuts(data::SBRPData, model, x, y)
  # setup
  B, A, n_cuts = data.B, data.D.A, 0
  # max clique model
  max_clique = direct_model(CPLEX.Optimizer())
  set_silent(max_clique)
  @variable(max_clique, z[block in B], Bin)
  @objective(max_clique, Max, ∑(z[block] for block in B))
  @constraint(max_clique, [(block, block′) in [(B[i], B[j]) for i in 1:length(B) for j in i + 1:length(B) if isempty(∩(B[i], B[j]))]], 1 >= z[block] + z[block′])
  @constraint(max_clique, ∑(z[block] for block in B) >= 2)
  while true
    # solve
    optimize!(max_clique)
    # base case
    termination_status(max_clique) == MOI.INFEASIBLE && break 
    # get clique
    B′ = filter(block -> value(z[block]) > 0.5, B)
    intersection = ∩(B′...)
    # check if it is a clique 
    println("Maximal clique: ", length(B′))
    println(B′)
    flush(stdout)
    all([!isempty(∩(block, block′)) for block′ ∈ B for block ∈ B′ if block′ ≠ block]) && error("It is not a clique")
    # add inquelity in the SBRP model
    @constraint(model, [block in B′], 1 - ∑(x[a] for a ∈ δ⁺(A, intersection)) >= ∑(∑(x[a] for a ∈ δ⁺(A, i)) for i ∈ setdiff(block, intersection, ∪(i for block′ ∈ B′ for i in ∩(block′, block) if block′ ∉ B′))))
    @constraint(model, [block in B′], ∑(∑(x[a] for a ∈ δ⁺(A, i) if a[2] ∈ intersection) for i ∈ setdiff(intersection, ∪(i for block′ ∈ B′ for i in ∩(block′, block) if block′ ∉ B′))) == 0)
    # update clique model
    @constraint(max_clique, ∑(z[block] for block ∈ B′) ≤ length(B′) - 1)
    # increase cuts number
    n_cuts += 1
  end
  return n_cuts
end

include("model_sbrp_max.jl")

include("model_sbrp_max_complete.jl")

end
