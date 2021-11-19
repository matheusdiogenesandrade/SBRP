module ModelSBRPComplete

using ..Model
using ...Data
using ...Data.SBRP
using CPLEX
using JuMP

export build_model_sbrp_complete

function build_model_sbrp_complete(data::SBRPData, app::Dict{String,Any})
  B, A, T, depot, Vb, gh_cuts, bfs_cuts = data.B, data.D.A, data.T, data.depot, Set{Int64}([i for b in data.B for i in b]), [], []
  V = Set{Int64}(vcat([i for (i, j) in A], [j for (i, j) in A]))
  # Formulation
  model = direct_model(CPLEX.Optimizer())
  set_silent(model)
#  MOI.set(model, MOI.NumberOfThreads(), 1)
  @variable(model, x[a in A], Bin)
  @variable(model, y[a in A], lower_bound=0, upper_bound=length(V))
  @objective(model, Min, sum(Data.time(data, a) * x[a] for a in A))
  @constraint(model, degree[i in V], sum(x[a] for a in δ⁻(A, i)) == sum(x[a] for a in δ⁺(A, i)))
  @constraint(model, sum(x[a] for a in δ⁺(A, depot)) <= 1)
  @constraint(model, mtz[i in Vb], sum(y[a] for a in δ⁺(A, i)) == sum(y[a] for a in δ⁻(A, i)) + sum(x[a] for a in δ⁺(A, i)))
  @constraint(model, mtz1[a in A], y[a] >= x[a])
  @constraint(model, mtz2[a in A], x[a] * length(V) >= y[a])
  @constraint(model, block[block in B], sum(sum(x[a] for a in δ⁺(A, i)) for i in block) == 1)
#  @constraint(model, sum(x[a] * time(data, a) for a in A) <= T - sum(time_block(data, block) for block in B))
  return (model, x, y)
end

function build_atsp_instance(data::SBRPData)
  # polynomial reduction 
  # SBRP attrs
  B, A, T, depot = data.B, data.D.A, data.T, data.depot
  V = Set{Int64}(vcat([i for (i, j) in A], [j for (i, j) in A]))
  # TSP attrs
  Vb, Vb′, n = Dict{Tuple{Int64, Array{Int64, 1}}, Int64}(), Dict{Tuple{Int64, Array{Int64, 1}}, Int64}(), max(collect(V)...) + 1
  # dummy nodes
  [(Vb[(i, b)] = n; n = n + 1) for b in B for i in b]
  [(Vb′[(i, b)] = n; n = n + 1) for b in B for i in b]
  # weights
  costs = merge!(
    Dict{Tuple{Int64, Int64}, Float64}((Vb[(i, b)],       Vb′[(i, b)])       => 0                       for b in B for i in b), # cycle arcs
    Dict{Tuple{Int64, Int64}, Float64}((Vb′[(b[i], b)],   Vb[(b[i + 1], b)]) => 0                       for b in B for i in 1:(length(b) - 1)), # cycle arcs
    Dict{Tuple{Int64, Int64}, Float64}((Vb′[(b[end], b)], Vb[(b[begin], b)]) => 0                       for b in B), # cycle arcs
    Dict{Tuple{Int64, Int64}, Float64}((Vb[(i, b)],       depot)             => 0                       for b in B for i in b), # depor arcs 
    Dict{Tuple{Int64, Int64}, Float64}((depot,            Vb′[(i, b)])       => 0                       for b in B for i in b), # depor arcs 
    Dict{Tuple{Int64, Int64}, Float64}((Vb[(i, b)],       Vb′[(j, b′)])      => Data.time(data, (i, j)) for b in B for b′ in B for i in b for j in b′ if b != b′ && (i, j) in keys(data.D.distance)) # block arcs
  )
  # mip atsp
  Aᵗ = collect(keys(costs))
  Vᵗ = collect(Set{Int64}(vcat([i for (i, j) in Aᵗ], [j for (i, j) in Aᵗ])))
  return Vᵗ, Aᵗ, costs, Vb, Vb′
end


end
