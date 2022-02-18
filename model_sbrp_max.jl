module ModelSBRPMax

include("symbols.jl")

using ..Model
using ...Data
using ...Data.SBRP
using CPLEX
using JuMP

export build_model_sbrp_max

include("SparseMaxFlowMinCut.jl")

add_subtour_cuts(model, sets::Set{Tuple{Arcs, Int}}) = add_cuts(model, [(∑(model[:x][a] for a in Aₛ), model[:z][i]) for (Aₛ, i) in sets])

function lazy_separation(data::SBRPData, Vb::Si, info, model, cb_data::CPLEX.CallbackContext, context_id::Clong)
  x, z = model[:x], model[:z]

  # preliminary checkings
  context_id != CPX_CALLBACKCONTEXT_CANDIDATE && return
  ispoint_p = Ref{Cint}()
  (CPXcallbackcandidateispoint(cb_data, ispoint_p) != 0 || ispoint_p[] == 0) && return 

  # get values
  CPLEX.load_callback_variable_primal(cb_data, context_id)
  x_val, z_val = callback_value.(Ref(cb_data), x), callback_value.(Ref(cb_data), z)

  # bfs
  A = data.D.A
  A′ = [a for a in A if x_val[a] > 0.5]
  sets = bfs_sets(A′, Vb, data.depot)

  # get valid components
#  components = Set{Tuple{Si, Int, Vi}}((S, i, block) for S in sets for block in blocks(data.B, S) for i in intersect(block, S) if z_val[i] + y_val[block] - 1 > length(δ⁺(A′, S)) + 0.5) 
  components = Set{Tuple{Arcs, Int}}()
  for S in sets 
    for i in S 
      if i in Vb && ∑(x_val[a] for a in δ⁺(A′, S)) + 0.5 < z_val[i]
        Aₛ = δ⁺(A, S)
        push!(components, (Aₛ, i))
        MOI.submit(model, MOI.LazyConstraint(cb_data), @build_constraint(Σ(x[a] for a in Aₛ) >= z[i])) 
      end
    end
  end

  # add ineqs
#  [MOI.submit(model, MOI.LazyConstraint(cb_data), @build_constraint(Σ(x[a] for a in δ⁺(A, S)) >= z[i] + y[block] - 1)) for (S, i, block) in components]

  # update info
  info["lazyCuts"] = info["lazyCuts"] + length(components)

  #
#  [flush_println(collect(S), i, block) for (S, i, block) in components]
end

function get_subtour_cuts(data::SBRPData, model, info::Dict{String, Any})

  # setup
  x, z = model[:x], model[:z]
  A, B, depot, components, Vb = data.D.A, data.B, data.depot, Set{Tuple{Arcs, Int}}(), blocks_nodes(data.B) 
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
    x_val = ArcCostMap(a => value(x[a]) for a in A)
    z_val = Dict{Int, Float64}(i => value(z[i]) for i in Vb)

    # get subsets
    g, M, newComponents = SparseMaxFlowMinCut.ArcFlow[], 100000, Set{Tuple{Arcs, Int}}()
    A′ = [a for a in A if x_val[a] > EPS]

    # mounting graph
    [push!(g, SparseMaxFlowMinCut.ArcFlow(Vₘ[a[1]], Vₘ[a[2]], trunc(floor(x_val[a], digits=5) * M))) for a in A′]

    # get subsets
    for source in Vb 

      # init
      maxFlow, flows, set = SparseMaxFlowMinCut.find_maxflow_mincut(SparseMaxFlowMinCut.Graph(n, g), Vₘ[source], Vₘ[depot])

      # base case 1
      (maxFlow / M >= z_val[source] - 1e-2 || set[Vₘ[depot]] == 1) && continue

      # get set
      S = Si([Vₘʳ[i] for i in 1:n if set[i] == 1])

      # base case 3
#      length(S) <= 1 && continue
      flush_println(length(S))

      # get components
#      [push!(newComponents, (S, i, block)) for block in blocks(B, S) for i in intersect(block, S) if z_val[i] + y_val[block] - 1 > flow + 1e-2]
      Aₛ = δ⁺(A, S)
      push!(newComponents, (Aₛ, source))
      add_cuts(model, [(∑(model[:x][a] for a in Aₛ), model[:z][source])])

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

subtour_cuts = Set{Tuple{Arcs, Int}}()

function build_model_sbrp_max(data::SBRPData, app::Dict{String,Any})
  B, A, T, V, profits, Vb, info = data.B, data.D.A, data.T, Set{Int}(vcat([i for (i, j) in data.D.A], [j for (i, j) in data.D.A])), data.profits, Set{Int}([i for b in data.B for i in b]), Dict{String, Any}("lazyCuts" => 0)

  function create_model(relax_x::Bool = false, relax_y::Bool = false, relax_z::Bool = false)
    global subtour_cuts
    model = direct_model(CPLEX.Optimizer())
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

    if relax_z
      @variable(model, z[i in Vb], lower_bound = 0, upper_bound = 1)
    else
      @variable(model, z[i in Vb], Bin)
    end

    @objective(model, Max, sum(data.profits[b] * y[b] for b in B))
    @constraint(model, degree[i in V], sum(x[a] for a in δ⁻(A, i)) == sum(x[a] for a in δ⁺(A, i)))
    @constraint(model, sum(x[a] for a in δ⁺(A, data.depot)) == 1)
    @constraint(model, block[block in B], sum(x[a] for i in block for a in δ⁺(A, i)) >= y[block])
    @constraint(model, sum(Data.SBRP.time(data, a) * x[a] for a in A) <= T - sum(y[block] * time_block(data, block) for block in B))
    @constraint(model, update_z[i in Vb], sum(x[a] for a in δ⁺(A, i)) <= (length(Vb)^2) * z[i])
    # subtour cuts
    add_subtour_cuts(model, subtour_cuts) 
    # registries
    model[:x], model[:y], model[:z] = x, y, z
    return model
  end

  # getting initial relaxation with all variables relaxed
  flush_println("Getting initial relaxation")
  model = create_model(true, true, true)
  optimize!(model)
  info["initialLP"] = objective_value(model)

  # getting initial relaxation with x and z relaxed (y integer)
  flush_println("Getting initial relaxation with y as integer")
  model = create_model(true, false, true)
  info["yLPTime"] = @elapsed optimize!(model)
  info["yLP"] = objective_value(model)

  # get max-flow cuts with x and y relaxed or integer
  model = create_model(true, !app["y-integer"])
  optimize!(model)
  info["maxFlowCutsTime"] = @elapsed subtour_cuts = get_subtour_cuts(data, model, info)
  info["maxFlowCuts"], info["maxFlowLP"] = length(subtour_cuts), objective_value(model)

  # integer model
  model = create_model(); 
  MOI.set(model, MOI.NumberOfThreads(), 1); # thread numbers
  MOI.set(model, CPLEX.CallbackFunction(), (cb_data, context_id) -> lazy_separation(data, Vb, info, model, cb_data, context_id)) # lazy callback

  # return
  return (model, model[:x], model[:y], model[:z], info)
end

end
