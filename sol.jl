using CPLEX
using DataStructures

function writesol(ids::Dict{Int64, Int64}, path::String, data::SBRPData, x, model, app::Dict{String, Any})
  tour = gettour(data, x)
  open(path, "w") do file
    write(file, "$(objective_value(model))\n")
    [write(file, "$node, ") for node in tour]
  end
end

function gettour(data::SBRPData, x)
  V, A, depot, tour = 1:length(data.D.V), [a for a in data.D.A if value(x[a]) > 0.5], data.depot, []
  #hierholzer's
  adjList, curr_path = Dict{Int64, Vector{Int64}}(i => [j for (i, j) in δ⁺(A, i) for i in 1:Int(floor(value(x[(i, j)]) + 0.5))] for i in V), Stack{Int}()
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
