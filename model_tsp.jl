module ModelTSP

using CPLEX
using JuMP
using ..Model
using ...Data
using ...Data.SBRP
#using CPLEX
#using JuMP

export solve_tsp, from_atsp_to_tsp, build_atsp_instance, build_model_atsp

#=
# The SBRPData must be a complete digraph where V = Vb \cup {0}
=#
function build_atsp_instance(data::SBRPData)
  # polynomial reduction 
  # SBRP attrs
  B, A, depot = data.B, data.D.A, data.depot
  V = Set{Int64}(vcat([i for (i, j) in A], [j for (i, j) in A]))
  # TSP attrs
  Vb, Vb′, n = Dict{Tuple{Int64, Array{Int64, 1}}, Int64}(), Dict{Tuple{Int64, Array{Int64, 1}}, Int64}(), max(collect(V)...) + 1
  # dummy nodes
  [(Vb[(i, b)] = n; n = n + 1) for b in B for i in b]
  [(Vb′[(i, b)] = n; n = n + 1) for b in B for i in b]
  # weights
  costs = merge!(
    Dict{Tuple{Int64, Int64}, Int64}((Vb[(i, b)],       Vb′[(i, b)])       => 0                                     for b in B for i in b), # cycle arcs
    Dict{Tuple{Int64, Int64}, Int64}((Vb′[(b[i], b)],   Vb[(b[i + 1], b)]) => 0                                     for b in B for i in 1:(length(b) - 1)), # cycle arcs
    Dict{Tuple{Int64, Int64}, Int64}((Vb′[(b[end], b)], Vb[(b[begin], b)]) => 0                                     for b in B), # cycle arcs
    Dict{Tuple{Int64, Int64}, Int64}((Vb[(i, b)],       depot)             => 0                                     for b in B for i in b), # depot arcs 
    Dict{Tuple{Int64, Int64}, Int64}((depot,            Vb′[(i, b)])       => 0                                     for b in B for i in b), # depot arcs 
    Dict{Tuple{Int64, Int64}, Int64}((Vb[(i, b)],       Vb′[(j, b′)])      => floor(Int64, data.D.distance[(i, j)]) for b in B for b′ in B for i in b for j in b′ if b != b′ && (i, j) in keys(data.D.distance)) # block arcs
  )
  # atsp
  Aᵗ = collect(keys(costs))
  Vᵗ = collect(Set{Int64}(vcat([i for (i, j) in Aᵗ], [j for (i, j) in Aᵗ])))
  return Vᵗ, Aᵗ, costs, Vb, Vb′
end

function from_atsp_to_tsp(V::Array{Int64, 1}, A::Array{Tuple{Int64, Int64}, 1}, costs::Dict{Tuple{Int64, Int64}, Int64})
  # convert to symmetric https://doi.org/10.1016/0167-6377(83)90048-2
  n = max(V...)
  V′, max_weight = vcat(V, [n + i for i in V]), 99999
  get_cost(i, j) = (i, j) in keys(costs) ? costs[(i, j)] : max_weight
  return [[(i <= n && j <= n) || (i > n && j > n) ? max_weight : (i > n && j <= n ? get_cost(i - n, j) : get_cost(j - n, i)) for i in V′] for j in V′]
end

function build_model_atsp(V, A, costs)
  # Formulation
  model = direct_model(CPLEX.Optimizer())
  set_silent(model)
  V′ = copy(V)
  pop!(V′)
  @variable(model, x[a in A], Bin)
  @variable(model, y[a in A], lower_bound=0, upper_bound=length(V))
  @objective(model, Min, sum(costs[a] * x[a] for a in A))
  @constraint(model, degree[i in V], sum(x[a] for a in δ⁻(A, i)) == sum(x[a] for a in δ⁺(A, i)))
  @constraint(model, once[i in V], sum(x[a] for a in δ⁻(A, i)) == 1)
  @constraint(model, mtz[i in V′], sum(y[a] for a in δ⁺(A, i)) == sum(y[a] for a in δ⁻(A, i)) + sum(x[a] for a in δ⁺(A, i)))
  @constraint(model, mtz1[a in A], y[a] >= x[a])
  @constraint(model, mtz2[a in A], x[a] * length(V) >= y[a])
  # return
  return (model, x, y)
end

function solve_tsp(costs::Array{Array{Int64, 1}, 1})
  concorde_path = "/home/matheusdiogenesandrade/concorde/TSP/concorde"
  n = length(costs)
  # write concorde input file
  name = "temp"
  filename = name * ".tsp"
  open(filename, "w") do io
    println(io, "NAME: $name")
    println(io, "TYPE: TSP")
    println(io, "COMMENT: $name")
    println(io, "DIMENSION:  $n")
    println(io, "EDGE_WEIGHT_TYPE: EXPLICIT")
    println(io, "EDGE_WEIGHT_FORMAT: FULL_MATRIX")
    println(io, "EDGE_WEIGHT_SECTION")
    for i in 1:length(costs)
      for j in 1:length(costs[i])
        print(io, "$(costs[i][j]) ")
      end
      println(io, "")
    end
    println(io, "EOF")
  end
  # run concorde
#  println(read(`$concorde_path $filename`, String))
  status = run(`$concorde_path $filename`, wait = false)
  while !success(status)
  end   
  # read solution
  tour = []
  open("$name.sol", "r") do sol
    #ignore first line
    lines = readlines(sol)
    popfirst!(lines)
    for line in lines
      parts = split(line, " ", keepempty = false)
      [push!(tour, parse(Int64, part) + 1) for part in parts] 
    end
  end
  # calculate cost
#  println(tour)
  cost = sum(costs[tour[i]][tour[i - 1]] for i in 2:length(tour)) + costs[tour[end]][tour[1]]
  # clean
  rm("$filename")
  rm("$name.sol")
  isfile("$name.mas") && rm("$name.mas")
  isfile("$name.sav") && rm("$name.sav")
  isfile("$name.pul") && rm("$name.pul")
  return cost
end

end
