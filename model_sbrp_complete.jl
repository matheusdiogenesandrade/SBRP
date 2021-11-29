module ModelSBRPComplete

using ..Model
using ...Data
using ...Data.SBRP
using CPLEX
using JuMP

export build_model_sbrp_complete

function build_model_sbrp_complete(data::SBRPData, app::Dict{String,Any})
  B, A, T, depot, Vb, info = data.B, data.D.A, data.T, data.depot, Set{Int64}([i for b in data.B for i in b]), Dict{String, Any}("lazyCuts" => 0)
  V = Set{Int64}(vcat([i for (i, j) in A], [j for (i, j) in A]))
  function add_basic_constraints(model)
    @objective(model, Min, sum(data.D.distance[a] * x[a] for a in A))
    @constraint(model, degree[i in V], sum(x[a] for a in δ⁻(A, i)) == sum(x[a] for a in δ⁺(A, i)))
    @constraint(model, sum(x[a] for a in δ⁺(A, depot)) <= 1)
    @constraint(model, mtz[i in Vb], sum(y[a] for a in δ⁺(A, i)) == sum(y[a] for a in δ⁻(A, i)) + sum(x[a] for a in δ⁺(A, i)))
    @constraint(model, mtz1[a in A], y[a] >= x[a])
    @constraint(model, mtz2[a in A], x[a] * length(V) >= y[a])
    @constraint(model, block[block in B], sum(x[a] for a in δ⁺(A, block)) == 1)
#    @constraint(model, sum(x[a] * time(data, a) for a in A) <= T - sum(time_block(data, block) for block in B))
  end
  # Formulation
  # frac model - get max-flow/min-cut cuts
  model = direct_model(CPLEX.Optimizer())
  set_silent(model)
  @variable(model, x[a in A], lower_bound = 0, upper_bound = 1)
  @variable(model, y[a in A], lower_bound=0, upper_bound=length(V))
  add_basic_constraints(model)
  info["maxFlowCutsTime"] = @elapsed sets = get_max_flow_min_cut_cuts(data, model, x)
  info["maxFlowCuts"], info["rootLP"] = length(sets), objective_value(model)
  # integer model
  model = direct_model(CPLEX.Optimizer())
  set_silent(model)
  MOI.set(model, MOI.NumberOfThreads(), 1)
  @variable(model, x[a in A], Bin)
  @variable(model, y[a in A], lower_bound=0, upper_bound=length(V))
  add_basic_constraints(model)
  add_subtour_ineqs(model, x, sets, A)
  # return
  return (model, x, y, info)
end

end
