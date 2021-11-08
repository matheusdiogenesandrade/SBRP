using CPLEX
using DataStructures

function writesol(ids::Dict{Int64, Int64}, path::String, data::SBRPData, x, model, app::Dict{String, Any})
  V, A, depot, tour = data.D.V, data.D.A, data.depot, []
  open(path, "w") do file
    write(file, "$(objective_value(model))\n")
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
    [write(file, "$node, ") for node in reverse(tour)]
  end
end
