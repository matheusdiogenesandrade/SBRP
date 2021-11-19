module ModelSBRP

using ..Model
using ...Data
using ...Data.SBRP
using CPLEX
using JuMP

export build_model_sbrp

function build_model_sbrp(data::SBRPData, app::Dict{String,Any})
  B, A, T, V, depot, Vb, bfs_cuts = data.B, data.D.A, data.T, 1:length(data.D.V), data.depot, Set{Int64}([i for b in data.B for i in b]), [], Set{Set{Int64}}()
  function add_basic_constraints(model)
    @objective(model, Min, sum(Data.time(data, a) * x[a] for a in A))
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

end
