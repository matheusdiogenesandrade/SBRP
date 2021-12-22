module ModelSBRPMaxComplete

using ..Model
using ...Data
using ...Data.SBRP
using CPLEX
using JuMP

export build_model_sbrp_max_complete

function build_model_sbrp_max_complete(data::SBRPData, app::Dict{String,Any})
  B, A, T, V, profits, Vb, info = data.B, data.D.A, data.T, Set{Int64}(vcat([i for (i, j) in data.D.A], [j for (i, j) in data.D.A])), data.profits, Set{Int64}([i for b in data.B for i in b]), Dict{String, Any}("lazyCuts" => 0)
  nodes_blocks = Dict{Int64, Array{Array{Int64}}}(i => [block for block in B if i in B] for i in Vb)
  function create_model(relax_x::Bool = false)
    model = direct_model(CPLEX.Optimizer())
    set_silent(model)
    if relax_x
      @variable(model, x[a in A], lower_bound = 0, upper_bound = 1)
    else
      @variable(model, x[a in A], Bin)
    end
    @variable(model, y[b in B], Bin)
    @objective(model, Max, sum(data.profits[b] * y[b] for b in B))
    @constraint(model, degree[i in V], sum(x[a] for a in δ⁻(A, i)) == sum(x[a] for a in δ⁺(A, i)))
    @constraint(model, degree_at_most_once[i in V], sum(x[a] for a in δ⁺(A, i)) <= 1)
    @constraint(model, block1[block in B], sum(x[a] for a in δ⁺(A, block)) >= y[block])
    @constraint(model, sum(Data.SBRP.time(data, a) * x[a] for a in A) <= T - sum(y[block] * time_block(data, block) for block in B))
    # improvements
#    @constraint(model, block2[block in B], sum(x[a] for a in δ⁺(A, block)) <= 1)
@constraint(model, block3[(i, block) in [(i, first(blocks)) for (i, blocks) in nodes_blocks if length(blocks) == 1]], y[block] - ∑(x[a] in vcat(δ⁺(A, i), δ⁻(A, i))) >= 0)
    @constraint(model, block4[block in [block for blocks in B if all([length(nodes_blocks[i]) == 1 for i in blocks])]], sum(x[(i, j)] for (i, j) in A if i in block && j in block) == 0)
    return model, x, y
  end
  # frac model - get max-flow/min-cut cuts and CC cuts
  model, x, y = create_model(true)
  # get initial LP
  optimize!(model); info["initialLP"] = objective_value(model)
  # add cuts
  info["maxFlowCutsTime"] = @elapsed sets_max_flow = get_max_flow_min_cut_cuts(data, model, x, y, info)
  #info["CCcutsTime"] = @elapsed sets_cc = get_cc_cuts(data, model, x, y, info)
  info["intersectionCutsTime"] = @elapsed info["intersectionCuts"] = add_intersection_cuts(data, model, x, y)
  # integer model
  model, x, y = create_model()
  MOI.set(model, MOI.NumberOfThreads(), 1)
  add_subtour_ineqs(model, x, y, sets_max_flow, A)
#  add_cc_ineqs(model, x, A, sets_cc)
  # lazy callback
  MOI.set(model, CPLEX.CallbackFunction(), (cb_data, context_id) -> lazy_separation(data, Vb, info, model, x, y, cb_data, context_id))
  return (model, x, y, info)
end

end
