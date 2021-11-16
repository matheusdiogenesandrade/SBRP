using CPLEX
using DataStructures

function writesol(ids::Dict{Int64, Int64}, path::String, data::SBRPData, x, model, app::Dict{String, Any})
  tour = gettour(data, x)
  open(path, "w") do file
    write(file, "$(objective_value(model))\n")
    [write(file, "$node, ") for node in tour]
  end
end

#function gettour(data::SBRPData, x)
function gettour(V::Array{Int64}, A::Array{Tuple{Int64, Int64}}, depot::Int64, x)
#  V, A, depot, tour = V, [a for a in data.D.A if value(x[a]) > 0.5], data.depot, []
  A′, tour = [a for a in A if value(x[a]) > 0.5], Array{Int64, 1}()
  #hierholzer's
  adjList, curr_path = Dict{Int64, Vector{Int64}}(i => [j for (i, j) in δ⁺(A′, i) for i in 1:Int(floor(value(x[(i, j)]) + 0.5))] for i in V), Stack{Int}()
  push!(curr_path, depot)
  curr_v = first(curr_path)
  while !isempty(curr_path)
    if !isempty(adjList[curr_v])
      push!(curr_path, curr_v)
      curr_v = pop!(adjList[curr_v])
    else
      push!(tour, curr_v)
      curr_v = pop!(curr_path)
    end
  end
  return reverse(tour)
end

function check_atsp_sol(tour::Array{Int64, 1}, Vb::Dict{Tuple{Int64, Array{Int64, 1}}, Int64}, Vb′::Dict{Tuple{Int64, Array{Int64, 1}}, Int64})
  Vₘ = merge(
    Dict{Int64, Array{Int64, 1}}(Vb[(i, b)] => b for (i, b) in keys(Vb)),
    Dict{Int64, Array{Int64, 1}}(Vb′[(i, b)] => b for (i, b) in keys(Vb′))
  )
  i, n = 1, length(tour)
  while i <= n
    node = tour[i]
    if node in keys(Vₘ)
      b, k = Vₘ[node], i + 1
      while k <= n && tour[k] in keys(Vₘ) && Vₘ[tour[k]] == b
        k = k + 1
      end
      k - i + 1 < length(b) * 2 && return false
      i = k
    else
      i = i + 1
    end
  end
  return true
end
