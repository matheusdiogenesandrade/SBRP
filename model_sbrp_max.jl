module ModelSBRPMax

using ..Model
using ...Data
using ...Data.SBRP
using CPLEX
using JuMP

export build_model_sbrp_max

function build_model_sbrp_max(data::SBRPData, app::Dict{String,Any})
  B, A, T, V, profits, Vb, info = data.B, data.D.A, data.T, Set{Int}(vcat([i for (i, j) in data.D.A], [j for (i, j) in data.D.A])), data.profits, Set{Int}([i for b in data.B for i in b]), Dict{String, Any}("lazyCuts" => 0)
  function create_model(relax_x::Bool = false)
    model = direct_model(CPLEX.Optimizer())
    set_silent(model)
    if relax_x
      @variable(model, x[a in A], lower_bound = 0, upper_bound = 1)
    else
      @variable(model, x[a in A], Int, lower_bound = 0, upper_bound = length(B))
    end
    @variable(model, y[b in B], Bin)
    @objective(model, Max, sum(data.profits[b] * y[b] for b in B))
    @constraint(model, degree[i in V], sum(x[a] for a in δ⁻(A, i)) == sum(x[a] for a in δ⁺(A, i)))
    @constraint(model, block1[block in B], sum(x[a] for a in δ⁺(A, block)) >= y[block])
    @constraint(model, sum(x[a] for a in δ⁺(A, data.depot)) <= 1)
    @constraint(model, sum(Data.time(data, a) * x[a] for a in A) <= T - sum(y[block] * time_block(data, block) for block in B))
    return model, x, y
  end
  # frac model - get max-flow/min-cut cuts
  model, x, y = create_model(true)
  info["maxFlowCutsTime"] = @elapsed sets = get_max_flow_min_cut_cuts(data, model, x, y, info)
  info["maxFlowCuts"], info["maxFlowLP"] = length(sets), objective_value(model)
  # integer model
  model, x, y = create_model()
  MOI.set(model, MOI.NumberOfThreads(), 1)
  add_subtour_ineqs(model, x, y, sets, A)
  # lazy callback
  MOI.set(model, CPLEX.CallbackFunction(), (cb_data, context_id) -> lazy_separation(data, Vb, info, model, x, y, cb_data, context_id))
  return (model, x, y, info)
end

end
