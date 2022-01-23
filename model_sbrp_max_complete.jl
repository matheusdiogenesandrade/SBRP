module ModelSBRPMaxComplete

include("symbols.jl")

using ..Model
using ...Data
using ...Data.SBRP
using CPLEX
using JuMP

export build_model_sbrp_max_complete

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
