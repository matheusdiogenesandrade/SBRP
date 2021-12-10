module SBRP

using ..Data
#import ..InputDigraph

export time_block, SBRPData, compact, readSBRPDataCarlos, readSBRPDataMatheus, distance_block, tour_distance, tour_time

mutable struct SBRPData
  D::Data.InputDigraph
  depot::Int64
  B::Array{Array{Int64, 1}}
  T::Float64
  profits::Dict{Array{Int64, 1}, Float64}
end

distance_block(data::SBRPData, block::Array{Int64, 1})                         = sum(data.D.distance[(block[i - 1], block[i])] for i in 2:length(block)) + data.D.distance[(block[end], block[begin])]
time_block(data::SBRPData, block::Array{Int64, 1})                             = 4 * distance_block(data, block) / NORMAL_SPEED
tour_distance(data::SBRPData, tour::Array{Int64, 1})                           = sum(data.D.distance[(tour[i - 1], tour[i])] for i in 2:length(tour))
tour_time(data::SBRPData, tour::Array{Int64, 1}, B::Array{Array{Int64, 1}, 1}) = ((tour_distance(data, tour) - sum(distance_block(data, block) for block in B)) / NORMAL_SPEED) + sum(time_block(data, block) for block in B)


function compact(data::SBRPData, V′)
  data′, V, paths = SBRPData(
    Data.InputDigraph(
                 Array{Vertex, 1}(vcat(data.D.V[data.depot], [data.D.V[i] for i in V′])), 
                 Array{Tuple{Int64,Int64}, 1}(), 
                 Dict{Tuple{Int64, Int64}, Float64}()
                ), 
    data.depot, data.B, data.T, # 2 hours
    data.profits
  ), 1:length(data.D.V), Dict{Tuple{Int64, Int64}, Array{Int64}}()
  for i in V′
    # bfs
    distances, pred, q = [typemax(Float64) for i in V], [i for i in V], [i]
    distances[i] = 0.0
    while !isempty(q)
      curr = popfirst!(q)
      for (curr, next) in δ⁺(data.D.A, curr)
        if !in(next, [data.depot, i]) && (pred[next] == next || distances[next] > distances[curr] + data.D.distance[(curr, next)])
          distances[next] = distances[curr] + data.D.distance[(curr, next)]
          pred[next] = curr
          push!(q, next)
        end
      end
    end
    # add cost and distance
    for j in V′
      pred[j] == j && continue
      data′.D.distance[(i, j)], paths[(i, j)], curr = distances[j], [], j
      while pred[curr] != curr
        pushfirst!(paths[(i, j)], curr)
        curr = pred[curr]
      end
      pushfirst!(paths[(i, j)], i)
    end
  end
  # add arcs
  [(data′.D.distance[(data.depot, i)] = data′.D.distance[(i, data.depot)] = 0.0) for i in V′]
  data′.D.A = collect(keys(data′.D.distance))
  # remove begins and ends of each path
  for (a, path) in paths
    popfirst!(path)
    pop!(path)
  end
  return data′, paths
end

function readSBRPDataCarlos(app::Dict{String,Any})
  depot = 1
  data, blocks, ids, blocks_profits = SBRPData(
    Data.InputDigraph(
                 Array{Vertex, 1}([Vertex(-1, -1, -1)]), 
                 Array{Tuple{Int64,Int64}, 1}(), 
                 Dict{Tuple{Int64, Int64}, Float64}()
                ), 
    depot,
    Array{Array{Int64, 1}, 1}(),
    120.0, # 2 hours
    Dict{Array{Int64, 1}, Float64}()
  ), Dict{Int64, Array{Tuple{Int64, Int64}, 1}}(), Dict{Int64, Int64}(), Dict{Int64, Int64}()
#  Vₘ = Dict{Int64, Int64}()
  open(app["instance"]) do f
    parts = split(readline(f), [' ']; limit=0, keepempty=false)
    # get params
    nNodes, nArcs, nBlocks = parse(Int64, parts[1]), parse(Int64, parts[2]), parse(Int64, parts[3])
    # get nodes
    for i in 1:nNodes
      parts = split(readline(f), [' ']; limit=0, keepempty=false)
      id = parse(Int64, parts[2])
      ids[id] = i + 1
      push!(data.D.V, Vertex(id, parse(Float64, parts[4]), parse(Float64, parts[3])))
    end
    # get distances
    for i in 1:nArcs
      parts = split(readline(f), [' ']; limit=0, keepempty=false)
      a = (ids[parse(Int64, parts[2])], ids[parse(Int64, parts[3])])
      data.D.distance[a] = floor(Int64, parse(Float64, parts[4]))
      # get blocks
      parts[5] == "-1" && continue
      id_block = parse(Int64, parts[5])
      !haskey(blocks, id_block) && (blocks[id_block] = Array{Tuple{Int64, Int64}, 1}())
      !haskey(blocks_profits, id_block) && (blocks_profits[id_block] = 0)
      push!(blocks[id_block], a)
      blocks_profits[id_block] += parse(Int64, parts[6])
    end
  end
  n_blocks, k = 16, 1
  # get blocks
  for (block, arcs) in blocks
    # get cycle     
    cycle, curr, next = Array{Int64, 1}(), first(arcs)[1], first(arcs)[2]
    push!(cycle, curr)
    while next != first(cycle)
      push!(cycle, next)
      curr, next = next, first(δ⁺(arcs, next))[2]
    end
    push!(data.B, cycle)
    data.profits[cycle] = blocks_profits[block]

#    k >= n_blocks && break
    k = k + 1
  end
  # add arcs
  data.D.A = [keys(data.D.distance)...]
  # compact
  Vb = Set{Int64}([i for b in data.B for i in b])
  data′, paths = compact(data, Vb)
  # update arcs
  data.D.A = collect(Set{Tuple{Int64, Int64}}(vcat(
    [(path[i], path[i + 1]) for (a, path) in paths for i in 1:(length(path) - 1)], # min paths arcs
    [(a[1], path[begin]) for (a, path) in paths if !isempty(path)], # min paths arcs
    [(path[end], a[2]) for (a, path) in paths if !isempty(path)], # min paths arcs
    [a for (a, path) in paths if isempty(path)], # min paths arcs (edge case)
    [(b[i], b[i + 1]) for b in data.B for i in 1:(length(b) - 1)], # blocks arcs
    [(b[end], b[begin]) for b in data.B], # blocks arcs
    [(depot, i) for i in Vb], # depot arcs
    [(i, depot) for i in Vb] # depot arcs
   )))
  # dummy weights
  [data.D.distance[(depot, i)] = data.D.distance[(i, depot)] = 0.0 for i in Vb]
  # check feasibility
  !check_sbrp_complete_feasibility(data′, Vb) && error("The SBRP instance is not feasible")
  # return
  return data, Dict{Int64, Int64}(v => k for (k, v) in ids), data′, paths
end

function readSBRPDataMatheus(app::Dict{String,Any})
  depot = 1
  data, blocks, ids = SBRPData(
    Data.InputDigraph(
                 Array{Vertex, 1}([Vertex(-1, -1, -1)]), 
                 Array{Tuple{Int64,Int64}, 1}(), 
                 Dict{Tuple{Int64, Int64}, Float64}()
                ), 
    depot,
    Array{Array{Int64, 1}, 1}(),
    120.0, # 2 hours
    Dict{Array{Int64, 1}, Float64}()
  ), Dict{Int64, Array{Tuple{Int64, Int64}, 1}}(), Dict{Int64, Int64}()
  open(app["instance"]) do f
    # ignore the first two lines
    readline(f)
    readline(f)
    nNodes, nArcs, nBlocks = parse(Int64, split(readline(f), [' ']; limit=0, keepempty=false)[end]), parse(Int64, split(readline(f), [' ']; limit=0, keepempty=false)[end]) + parse(Int64, split(readline(f), [' ']; limit=0, keepempty=false)[end]), parse(Int64, split(readline(f), [' ']; limit=0, keepempty=false)[end])
    readline(f)
    readline(f)
    # get arcs
    for n in 1:nArcs
      parts = split(readline(f), ['(', ')', ' ', ',']; limit=0, keepempty=false)
      data.D.distance[(parse(Int64, parts[1]), parse(Int64, parts[2]))] = floor(Int64, parse(Float64, parts[end]))
    end
    # get vertices
    readline(f)
    for i in 1:nNodes
      parts = split(readline(f), [' ']; limit=0, keepempty=false)
      id = parse(Int64, parts[1])
      ids[id] = i + 1
      push!(data.D.V, Vertex(id, parse(Float64, parts[2]), parse(Float64, parts[3])))
    end
    # get blocks
    readline(f)
    n_blocks, k = 2, 1
    for i in 1:nBlocks
      parts = split(readline(f), [',', ' ']; limit=0, keepempty=false)
      push!(data.B, Array{Int64, 1}([ids[parse(Int64, part)] for part in parts]))
#      k >= n_blocks && break
      k = k + 1
    end
  end
  # update arcs ids
  data.D.distance = Dict{Tuple{Int64, Int64}, Float64}((ids[a[1]], ids[a[2]]) => distance for (a, distance) in data.D.distance)
  data.D.A = collect(keys(data.D.distance))
  # compact
  Vb = Set{Int64}([i for b in data.B for i in b])
  data′, paths = compact(data, Vb)
  # update arcs
  data.D.A = collect(Set{Tuple{Int64, Int64}}(vcat(
    [(path[i], path[i + 1]) for (a, path) in paths for i in 1:(length(path) - 1)], # min paths arcs
    [(a[1], path[begin]) for (a, path) in paths if !isempty(path)], # min paths arcs
    [(path[end], a[2]) for (a, path) in paths if !isempty(path)], # min paths arcs
    [a for (a, path) in paths if isempty(path)], # min paths arcs (edge case)
    [(b[i], b[i + 1]) for b in data.B for i in 1:(length(b) - 1)], # blocks arcs
    [(b[end], b[begin]) for b in data.B], # blocks arcs
    [(depot, i) for i in Vb], # depot arcs
    [(i, depot) for i in Vb] # depot arcs
   )))
  # dummy weights
  [data.D.distance[(depot, i)] = data.D.distance[(i, depot)] = 0.0 for i in Vb]
  # check feasibility
  !check_sbrp_complete_feasibility(data′, Vb) && error("The Complete SBRP instance is not feasible")
  # define profits
  [data.profits[b] = 1.0 for b in data.B]
  # return
  return data, Dict{Int64, Int64}(v => k for (k, v) in ids), data′, paths
end

function check_sbrp_complete_feasibility(data_complete::SBRPData, Vb)
  A′ = keys(data_complete.D.distance)
  return all((i, j) in A′ && (j, i) in A′ for i in Vb for j in Vb if i < j)
end

end
