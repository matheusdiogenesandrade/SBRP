module ModelSBRPMaxComplete

include("symbols.jl")

using ..Model
using ...Data
using ...Data.SBRP
using CPLEX
using JuMP

export build_model_sbrp_max_complete

include("SparseMaxFlowMinCut.jl")

function add_intersection_cuts1(model, cuts::Vector{Arcs})
  add_cuts(model, [(sum(model[:x][a] for a in arcs), 0) for arcs in cuts])
  add_cuts(model, [(0, sum(model[:x][a] for a in arcs)) for arcs in cuts])
end

add_intersection_cuts2(model, cuts::Vector{Arcs}) = add_cuts(model, [(2, sum(model[:x][a] for a in arcs)) for arcs in cuts])

add_subtour_cuts(model, sets::Set{Tuple{Arcs, Arcs}}) = add_cuts(model, [(∑(model[:x][a] for a in Aₛ), ∑(∑(model[:x][a] for a in Aᵢ))) for (Aₛ, Aᵢ) in sets])

function lazy_separation(data::SBRPData, Vb, info, model, cb_data::CPLEX.CallbackContext, context_id::Clong)
  x, y = model[:x], model[:y]

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
#  components = Set{Tuple{Si, Int, Vi}}((S, i, block) for S in sets for block in blocks(data.B, S) for i in intersect(block, S) if z_val[i] + y_val[block] - 1 > length(δ⁺(A′, S)) + 0.5) 
  components = Set{Tuple{Arcs, Arcs}}()
  for S in sets 
    lhs = ∑(x_val[a] for a in δ⁺(A′, S))
    Aₛ, A′ₛ = δ⁺(A, S), δ⁺(A′, S) 
    for i in S 
      A′ᵢ = δ⁺(A′, i)
      if i in Vb && ∑(x_val[a] for a in A′ₛ) + 0.5 < ∑(x_val[a] for a in A′ᵢ)
        Aᵢ = δ⁺(A, i)
        push!(components, (Aₛ, Aᵢ))
        MOI.submit(model, MOI.LazyConstraint(cb_data), @build_constraint(Σ(x[a] for a in Aₛ) >= Σ(x[a] for a in Aᵢ))) 
      end
    end
  end

  # add ineqs
#  [MOI.submit(model, MOI.LazyConstraint(cb_data), @build_constraint(Σ(x[a] for a in δ⁺(A, S)) >= z[i] + y[block] - 1)) for (S, i, block) in components]

  # update info
  info["lazyCuts"] = info["lazyCuts"] + length(components)
end

function get_subtour_cuts(data::SBRPData, model, info::Dict{String, Any})
  x, y = model[:x], model[:y]

  # setup
  A, B, depot, components, Vb = data.D.A, data.B, data.depot, Set{Tuple{Arcs, Arcs}}(), blocks_nodes(data.B) 
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
    g, M, newComponents = SparseMaxFlowMinCut.ArcFlow[], 100000, Set{Tuple{Arcs, Arcs}}()
    A′ = [a for a in A if x_val[a] > EPS]

    # mounting graph
    [push!(g, SparseMaxFlowMinCut.ArcFlow(Vₘ[a[1]], Vₘ[a[2]], trunc(floor(x_val[a], digits=5) * M))) for a in A′]

    # get subsets
    for source in Vb 

      # init
      maxFlow, flows, set = SparseMaxFlowMinCut.find_maxflow_mincut(SparseMaxFlowMinCut.Graph(n, g), Vₘ[source], Vₘ[depot])
      flow = maxFlow / M

      # base case 1
      (flow >= 1 - 1e-2 || set[Vₘ[depot]] == 1) && continue

      # base case 2
      A′ᵢ = δ⁺(A′, source)
      rhs = ∑(x_val[a] for a in A′ᵢ)
        # get set
      S = Si([Vₘʳ[i] for i in 1:n if set[i] == 1])

      flow + 1e2 >= rhs && continue

      # base case 3
#      length(S) <= 1 && continue
      flush_println(length(S))

      # get components
#      [push!(newComponents, (S, i, block)) for block in blocks(B, S) for i in intersect(block, S) if z_val[i] + y_val[block] - 1 > flow + 1e-2]
      Aₛ = δ⁺(A, S)
      push!(newComponents, (Aₛ, Aᵢ))
      add_cuts(model, [(∑(x[a] for a in Aₛ), Σ(x[a] for a in Aᵢ))])

    end

    # base case
    isempty(newComponents) && break
#    flush_println(length(sets′))
#    [flush_println(S, i, block) for (S, i, block) in sets′]

    # update infos
#    info["iteration_" * string(iteration) * "_time"], info["iteration_" * string(iteration) * "_cuts"], iteration = time, length(sets′), iteration + 1

    # store components
    union!(components, newComponents)

    # add ineqs
    add_subtour_cuts(model, components)
  end
  return components
end

function get_intersection_cuts(data::SBRPData)

  # setup
  B, A, cuts1, cuts2, cliques = data.B, data.D.A, Vector{Arcs}(), Vector{Arcs}(), Vector{VVi}()

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
    # store clique
    push!(cliques, B′)
    # update clique model
    @constraint(max_clique, ∑(z[block] for block ∈ B′) ≤ length(B′) - 1)
  end

  # get nodes of each block
  Vb = blocks_nodes(B)
  nodes_blocks = Dict{Int, VVi}(i => [block for block in B if i in block] for i in Vb)

  # get cuts
  for clique in cliques
    # get intersection
    intersection = ∩(clique...)

    ######## intersection cuts 1 ########
    
    isolated_intersection = setdiff(intersection, ∪(i for block ∈ setdiff(B, clique) for i in block))
#    isolated_intersection = sediff(intersection, ∪(block... for block ∈ setdiff(B, clique)))
    length(isolated_intersection) > 1 && push!(cuts1, Arcs([(i, j) for (i, j) in χ(isolated_intersection) if i ≠ j]))

    flush_println("Isolated clique: ", isolated_intersection)

    ######## intersection cuts 2 ########
    
    flush_println("Clique: ", length(clique))

    for block in clique
      # get arcs incident to nodes that belong exclusively to the block 
#      independent_arcs = [a for (i, blocks) in nodes_blocks for a in δ⁺(A, i) if ∧(length(blocks) == 1, block ∈ blocks)]
      covered_arcs = [a for i in block for a in δ⁺(A, i) if ⊆(nodes_blocks[i], clique)]

      # edge case
#      isempty(independent_arcs) && continue
      isempty(covered_arcs) && continue

      # store
#      push!(cuts2, Arcs(∪([a for i ∈ intersection for a ∈ δ⁺(A, i)], independent_arcs)))
      push!(cuts2, Arcs(covered_arcs))
    end

  end

  return cuts1, cuts2
end

intersection_cuts1, intersection_cuts2 = Vector{Arcs}(), Vector{Arcs}()
subtour_cuts = Set{Tuple{Arcs, Arcs}}()

function build_model_sbrp_max_complete(data::SBRPData, app::Dict{String,Any})
    global STARTING_TIME 
    global COST_PER_TIME

  depot, B, A, T, V, profits, Vb, info = data.depot, data.B, data.D.A, data.T, data.D.V, data.profits, blocks_nodes(data.B), Dict{String, Any}("lazyCuts" => 0)
  A′, V′ = filter(a -> a[2] != depot, A), setdiff(keys(V), depot)
  P′ = filter(a -> a[1] < a[2], χ(V′))

  nodes_blocks = Dict{Int, VVi}(i => [block for block in B if i in block] for i in Vb)

  function create_model(relax_x::Bool = false, relax_y::Bool = false)
    global intersection_cuts1, intersection_cuts2, subtour_cuts
    model = direct_model(CPLEX.Optimizer())
    set_parameters(model, "CPX_PARAM_TILIM" => 3600)
#    set_silent(model)
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
    @constraint(model, degree[i in keys(V)], sum(x[a] for a in δ⁻(A, i)) == sum(x[a] for a in δ⁺(A, i)))
    @constraint(model, degree_at_most_once[i in keys(V)], sum(x[a] for a in δ⁺(A, i)) <= 1)
    @constraint(model, sum(x[a] for a in δ⁺(A, depot)) == 1)
    @constraint(model, block1[block in B], sum(x[a] for a in δ⁺(A, block)) >= y[block])
    @constraint(model, sum(Data.SBRP.time(data, a) * x[a] for a in A) <= T - sum(y[block] * time_block(data, block) for block in B))
    # improvements
    @constraint(model, block3[block in B], y[block] - ∑(x[a] in δ⁺(A, i) for i in block if length(nodes_blocks) == 1) >= 0)
    @constraint(model, block4[block in B], sum(x[(i, j)] for (i, j) in A if i in block && j in block && (length(nodes_blocks[i]) == 1 || length(nodes_blocks[j]) == 1)) == 0)
    @constraint(model, subcycle_size_two[(i, j) in P′], x[(i, j)] + x[(j, i)] <= 1)
    # MTZs
    if app["arcs-mtz"] # Arcs MTZ
      @variable(model, t[a in A], lower_bound = 0, upper_bound = T)
      @constraint(model, sum(t[a] for a in δ⁺(A, depot)) == 0.0)
      @constraint(model, mtz[i in V′], sum(t[a] for a in δ⁺(A, i)) == sum(t[a] for a in δ⁻(A, i)) + sum(x[a] * Data.SBRP.time(data, a) for a in δ⁺(A, i)))
      @constraint(model, ub[i in V′], sum(t[a] for a in δ⁻(A, i)) <= T - sum(y[block] * time_block(data, block) for block in B))
      @constraint(model, ub1[i in V′], sum(t[a] for a in δ⁻(A, i)) <= sum(x[a] for a in δ⁻(A, i)) * T)
      @constraint(model, ub2[a in A], t[a] <= x[a] * T)
    else # Nodes MTZ
      @variable(model, t[i in keys(V)], lower_bound = 0, upper_bound = T)
      @constraint(model, t[depot] == 0.0)
      @constraint(model, mtz[(i, j) in A′], t[j] >= t[i] + x[(i, j)] * Data.SBRP.time(data, (i, j)) - (1 - x[(i, j)]) * T - x[(j, i)] * Data.SBRP.time(data, (j, i)))
      @constraint(model, ub1[i in V′], t[i] <= T - sum(y[block] * time_block(data, block) for block in B))
      @constraint(model, ub2[i in V′], t[i] <= sum(x[a] for a in δ⁺(A, i)) * T)
      @constraint(model, ub3[i in V′], t[i] <= sum(x[a] for a in δ⁻(A, i)) * T)
    end
    # subtour cuts
    add_subtour_cuts(model, subtour_cuts) 
    # registries
    model[:x], model[:y] = x, y
    return model, x, y
  end
  # creating model with both variables (x and y) relaxed
  model, x, y = create_model(true, true)

  # getting intersection cuts
  if app["intersection-cuts"]
    flush_println("Getting intersection cuts")
    info["intersectionCutsTime"] = @elapsed intersection_cuts1, intersection_cuts2 = get_intersection_cuts(data) # get intersection cuts
    info["intersectionCuts1"], info["intersectionCuts2"] = length(intersection_cuts1), length(intersection_cuts2)
    # adding intersection cuts to the recently create model
    add_intersection_cuts1(model, intersection_cuts1)
    add_intersection_cuts2(model, intersection_cuts2)
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
  if app["subcycle-separation"]
    model, x, y = create_model(true, !app["y-integer"])
    optimize!(model)
    info["maxFlowCutsTime"] = @elapsed subtour_cuts = get_subtour_cuts(data, model, info)
    info["maxFlowCuts"], info["maxFlowLP"] = length(subtour_cuts), objective_value(model)
  end

  # integer model
  model, x, y = create_model()
#  MOI.set(model, MOI.NumberOfThreads(), 1) # thread numbers
#  MOI.set(model, CPLEX.CallbackFunction(), (cb_data, context_id) -> lazy_separation(data, Vb, info, model, cb_data, context_id)) # lazy callback
    MOI.set(model, CPLEX.CallbackFunction(), (cb_data, context_id) -> updateCostPerTimeRelation(data, model, cb_data, context_id)) # lazy callback

  # return
  return (model, x, y, info)
end


end
