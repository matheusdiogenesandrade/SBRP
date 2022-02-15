module ModelSBRPMaxComplete

include("symbols.jl")

using ..Model
using ...Data
using ...Data.SBRP
using CPLEX
using JuMP

export build_model_sbrp_max_complete

include("SparseMaxFlowMinCut.jl")

add_intersection_cuts(model, cuts::Vector{Tuple{Arcs, Arcs}}) = add_cuts(model, [(- 1 + sum(model[:x][a] for a in lhs), - sum(model[:x][a] for a in rhs)) for (lhs, rhs) in cuts])

add_subtour_cuts(model, x, y, sets::Set{Tuple{Set{Int}, Int, Vi}}, A::Arcs) = add_cuts(model, [(sum(x[a] for a in δ⁺(A, S)), sum(x[(i, j)] for (i, j) in δ⁺(A, i) if j in S) + y[block] - 1) for (S, i, block) in sets])

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
  sets′ = Set{Tuple{Set{Int}, Int, Vi}}((S, i, block) for S in sets for block in blocks(data.B, S) for i in intersect(block, S) if y_val[block] + sum(x_val[(i, j)] for (i, j) in δ⁺(A, i) if j in S && x_val[(i, j)] > EPS) - 1 > length(δ⁺(A′, S)) + 0.5)

  # add ineqs
  [MOI.submit(model, MOI.LazyConstraint(cb_data), @build_constraint(sum(x[a] for a in δ⁺(A, S)) >= y[block] + sum(x[(i, j)] for (i, j) in δ⁺(A, i) if j in S) - 1)) for (S, i, block) in sets′]

  # update info
  info["lazyCuts"] = info["lazyCuts"] + length(sets′)

  #
#  [flush_println(collect(S), i, block) for (S, i, block) in sets′]
end

function get_subtour_cuts(data::SBRPData, model, x, y, info::Dict{String, Any})

  # setup
  A, B, depot, components, Vb = data.D.A, data.B, data.depot, Set{Tuple{Si, Int, Vi}}(), Si([i for b in data.B for i in b])
  n, iteration = max(Si(vcat([i for (i, j) in A], [j for (i, j) in A]))...), 1
  Vₘ = Dict{Int, Int}(i => idx for (idx, i) in enumerate(keys(data.D.V)))
  Vₘʳ = Dict{Int, Int}(idx => i for (idx, i) in enumerate(keys(data.D.V)))
  n, iteration = length(Vₘ), 1

  # get cuts greedly
  while true
    # store time
    time = @elapsed optimize!(model)

    # error checking
    !in(termination_status(model), [MOI.OPTIMAL, MOI.TIME_LIMIT, MOI.ALMOST_INFEASIBLE]) && error("The model could not be solved")

    # get values
    y_val = Dict{Vi, Float64}(block => value(y[block]) for block in B)
    x_val = ArcCostMap(a => value(x[a]) for a in A)

    # get subsets
    g, M, newComponents = SparseMaxFlowMinCut.ArcFlow[], 100000, Set{Tuple{Si, Int, Vi}}()

    # mouting graph
    [push!(g, SparseMaxFlowMinCut.ArcFlow(Vₘ[a[1]], Vₘ[a[2]], trunc(floor(x_val[a], digits=5) * M))) for a in A if x_val[a] > EPS]

    # get subsets
    for source in Vb 

      # init
      maxFlow, flows, set = SparseMaxFlowMinCut.find_maxflow_mincut(SparseMaxFlowMinCut.Graph(n, g), Vₘ[source], Vₘ[depot])
      flow = maxFlow / M

      # base case 1
      (flow >= 1 - 1e-2 || set[Vₘ[depot]] == 1) && continue

      # get set
      S = Si([Vₘʳ[i] for i in 1:n if set[i] == 1])

      # base case 2
      length(S) <= 1 && continue
#      println(length(S) == length(Vb))

      # get components
      [push!(newComponents, (S, i, block)) for block in blocks(B, S) for i in intersect(block, S) if y_val[block] + sum(x_val[(i, j)] for (i, j) in δ⁺(A, i) if j in S) - 1 > flow + 1e-2]
    end

    # base case
    isempty(newComponents) && break
#    flush_println(length(newComponents))
#    [flush_println(S, i, block) for (S, i, block) in newComponents]

    # update infos
#    info["iteration_" * string(iteration) * "_time"], info["iteration_" * string(iteration) * "_cuts"], iteration = time, length(sets′), iteration + 1

    # store sets″
    union!(components, newComponents)

    # add ineqs
    add_subtour_cuts(model, x, y, newComponents, A)
  end
  return components 
end

function get_intersection_cuts(data::SBRPData, model, x, y)
  # setup
  B, A, cuts = data.B, data.D.A, Vector{Tuple{Arcs, Arcs}}()
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
    # check if it is a clique 
    all([!isempty(∩(block, block′)) for block′ ∈ B for block ∈ B′ if block′ ≠ block]) && error("It is not a clique")
    # get intersection
    intersection = ∩(B′...)
    flush_println("Maximal clique: ", length(B′))
    # add inquelity in the SBRP model
    for block in B′
      rhs = Arcs([a for i ∈ setdiff(block, intersection, ∪(i for block′ ∈ B′ for i in ∩(block′, block) if block′ ∉ B′)) for a ∈ δ⁺(A, i)])
      isempty(rhs) && continue
      push!(cuts, (Arcs([a for a ∈ δ⁺(A, intersection)]), rhs))
      flush_println(block)
    end
    @constraint(model, [block in B′], ∑(∑(x[a] for a ∈ δ⁺(A, i) if a[2] ∈ intersection) for i ∈ setdiff(intersection, ∪(i for block′ ∈ B′ for i in ∩(block′, block) if block′ ∉ B′))) == 0)
    # update clique model
    @constraint(max_clique, ∑(z[block] for block ∈ B′) ≤ length(B′) - 1)
    # increase cuts number
  end
  return cuts
end

intersection_cuts = Vector{Tuple{Arcs, Arcs}}()
subtour_cuts = Set{Tuple{Set{Int}, Int, Vi}}()

function build_model_sbrp_max_complete(data::SBRPData, app::Dict{String,Any})
  B, A, T, V, profits, Vb, info = data.B, data.D.A, data.T, Set{Int}(vcat([i for (i, j) in data.D.A], [j for (i, j) in data.D.A])), data.profits, Set{Int}([i for b in data.B for i in b]), Dict{String, Any}("lazyCuts" => 0)

  nodes_blocks = Dict{Int, Vector{Vi}}(i => [block for block in B if i in B] for i in Vb)

  function create_model(relax_x::Bool = false, relax_y::Bool = false)
    global intersection_cuts, subtour_cuts
    model = direct_model(CPLEX.Optimizer())
    set_silent(model)
    if relax_x
      @variable(model, x[a in A], lower_bound = 0, upper_bound = 1)
    else
      @variable(model, x[a in A], Bin)
    end
    if relax_y
      @variable(model, y[b in B], lower_bound = 0, upper_bound = 1)
    else
      @variable(model, y[b in B], Bin)
    end
    @objective(model, Max, sum(data.profits[b] * y[b] for b in B))
    @constraint(model, degree[i in V], sum(x[a] for a in δ⁻(A, i)) == sum(x[a] for a in δ⁺(A, i)))
    @constraint(model, degree_at_most_once[i in V], sum(x[a] for a in δ⁺(A, i)) <= 1)
    @constraint(model, block1[block in B], sum(x[a] for a in δ⁺(A, block)) >= y[block])
    @constraint(model, sum(Data.SBRP.time(data, a) * x[a] for a in A) <= T - sum(y[block] * time_block(data, block) for block in B))
    # improvements
    @constraint(model, block3[block in B], y[block] - ∑(x[a] in δ⁺(A, i) for i in block if length(nodes_blocks) == 1) >= 0)
    @constraint(model, block4[block in B], sum(x[(i, j)] for (i, j) in A if i in block && j in block && nodes_blocks[i] == nodes_blocks[j] && length(nodes_blocks[j]) == 1) == 0)
    # intersection cuts
    add_intersection_cuts(model, intersection_cuts)
    # subtour cuts
    add_subtour_cuts(model, x, y, subtour_cuts, A) 
    # registries
    model[:x], model[:y] = x, y
    return model, x, y
  end
  # creating model with both variables (x and y) relaxed
  model, x, y = create_model(true, true)

  # getting intersection cuts
  if app["intersection-cuts"]
    flush_println("Getting intersection cuts")
    info["intersectionCutsTime"] = @elapsed intersection_cuts = get_intersection_cuts(data, model, x, y) # get intersection cuts
    info["intersectionCuts"] = length(intersection_cuts)
    # adding intersection cuts to the recently create model
    add_intersection_cuts(model, intersection_cuts)
  end

  # getting initial relaxation with both variables <x and y> relaxed
  flush_println("Getting initial relaxation")
  optimize!(model)
  info["initialLP"] = objective_value(model)

  # getting initial relaxation with only x relaxed (y integer)
  flush_println("Getting initial relaxation with y as integer")
  model, x, y = create_model(true)
  info["yLPTime"] = @elapsed optimize!(model)
  info["yLP"] = objective_value(model)

  # get max-flow cuts with x and y relaxed or integer
  model, x, y = create_model(true, !app["y-integer"])
  optimize!(model)
  info["maxFlowCutsTime"] = @elapsed subtour_cuts = get_subtour_cuts(data, model, x, y, info)
  info["maxFlowCuts"], info["maxFlowLP"] = length(subtour_cuts), objective_value(model)

  # integer model
  model, x, y = create_model(); 
  MOI.set(model, MOI.NumberOfThreads(), 1); # thread numbers
  MOI.set(model, CPLEX.CallbackFunction(), (cb_data, context_id) -> lazy_separation(data, Vb, info, model, x, y, cb_data, context_id)) # lazy callback

  # return
  return (model, x, y, info)
end

end
