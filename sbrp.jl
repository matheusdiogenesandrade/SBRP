module SBRP

include("symbols.jl")

using ..Data
#import ..InputDigraph

export time_block, SBRPData, compact, readSBRPDataCarlos, readSBRPDataMatheus, distance_block, tour_distance, tour_time

mutable struct SBRPData
  D::Data.InputDigraph
  depot::Int
  B::Vector{Vi}
  T::Float64
  profits::Dict{Vi, Float64}
end

# 40 km / 60 min = 40 km/h
NORMAL_SPEED = (40.0 * 1e3)/60.0

time(data, a) = a[1] == a[2] ? 0.0 : data.D.distance[a] / NORMAL_SPEED
distance_block(data::SBRPData, block::Vi)                         = sum(data.D.distance[(block[i - 1], block[i])] for i in 2:length(block)) + data.D.distance[(block[end], block[begin])]
time_block(data::SBRPData, block::Vi)                             = 4 * distance_block(data, block) / NORMAL_SPEED
tour_distance(data::SBRPData, tour::Vi)                           = sum(data.D.distance[(tour[i - 1], tour[i])] for i in 2:length(tour))
tour_time(data::SBRPData, tour::Vi, B::Vector{Vi}) = ((tour_distance(data, tour) - sum(distance_block(data, block) for block in B)) / NORMAL_SPEED) + sum(time_block(data, block) for block in B)


function compact(data::SBRPData, V′)
  data′, V, paths = SBRPData(
    Data.InputDigraph(
                 Vector{Vertex}(vcat(data.D.V[data.depot], [data.D.V[i] for i in V′])), 
                 Arcs(), 
                 Dict{Arc, Float64}()
                ), 
    data.depot, data.B, data.T, # 2 hours
    data.profits
  ), 1:length(data.D.V), Dict{Tuple{Int, Int}, Vi}()
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
                 Vector{Vertex, 1}([Vertex(-1, -1, -1)]), 
                 Vector{Tuple{Int,Int}, 1}(), 
                 Dict{Tuple{Int, Int}, Float64}()
                ), 
    depot,
    Vector{Vi, 1}(),
    120.0, # 2 hours
    Dict{Vi, Float64}()
  ), Dict{Int, Vector{Tuple{Int, Int}, 1}}(), Dict{Int, Int}(), Dict{Int, Int}()
#  Vₘ = Dict{Int, Int}()
  open(app["instance"]) do f
    parts = split(readline(f), [' ']; limit=0, keepempty=false)
    # get params
    nNodes, nArcs, nBlocks = parse(Int, parts[1]), parse(Int, parts[2]), parse(Int, parts[3])
    # get nodes
    for i in 1:nNodes
      parts = split(readline(f), [' ']; limit=0, keepempty=false)
      id = parse(Int, parts[2])
      ids[id] = i + 1
      push!(data.D.V, Vertex(id, parse(Float64, parts[4]), parse(Float64, parts[3])))
    end
    # get distances
    for i in 1:nArcs
      parts = split(readline(f), [' ']; limit=0, keepempty=false)
      a = (ids[parse(Int, parts[2])], ids[parse(Int, parts[3])])
      data.D.distance[a] = floor(Int, parse(Float64, parts[4]))
      # get blocks
      parts[5] == "-1" && continue
      id_block = parse(Int, parts[5])
      !haskey(blocks, id_block) && (blocks[id_block] = Vector{Tuple{Int, Int}, 1}())
      !haskey(blocks_profits, id_block) && (blocks_profits[id_block] = 0)
      push!(blocks[id_block], a)
      blocks_profits[id_block] += parse(Int, parts[6])
    end
  end
  n_blocks, k = 16, 1
  # get blocks
  for (block, arcs) in blocks
    # get cycle     
    cycle, curr, next = Vi(), first(arcs)[1], first(arcs)[2]
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
  Vb = Set{Int}([i for b in data.B for i in b])
  data′, paths = compact(data, Vb)
  # update arcs
  data.D.A = collect(Set{Tuple{Int, Int}}(vcat(
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
  return data, Dict{Int, Int}(v => k for (k, v) in ids), data′, paths
end

function readSBRPDataMatheus(app::Dict{String,Any})
  depot = 1
  data, blocks, ids = SBRPData(
    Data.InputDigraph(
                 Vector{Vertex}([Vertex(-1, -1, -1)]), 
                 Arcs(), 
                 Dict{Arc, Float64}()
                ), 
    depot,
    Vector{Vi}(),
    120.0, # 2 hours
    Dict{Vi, Float64}()
  ), Dict{Int, Arcs}(), Dict{Int, Int}()
  open(app["instance"]) do f
    # ignore the first two lines
    readline(f)
    readline(f)
    nNodes, nArcs, nBlocks = parse(Int, split(readline(f), [' ']; limit=0, keepempty=false)[end]), parse(Int, split(readline(f), [' ']; limit=0, keepempty=false)[end]) + parse(Int, split(readline(f), [' ']; limit=0, keepempty=false)[end]), parse(Int, split(readline(f), [' ']; limit=0, keepempty=false)[end])
    readline(f)
    readline(f)
    # get arcs
    for n in 1:nArcs
      parts = split(readline(f), ['(', ')', ' ', ',']; limit=0, keepempty=false)
      data.D.distance[(parse(Int, parts[1]), parse(Int, parts[2]))] = floor(Int, parse(Float64, parts[end]))
    end
    # get vertices
    readline(f)
    for i in 1:nNodes
      parts = split(readline(f), [' ']; limit=0, keepempty=false)
      id = parse(Int, parts[1])
      ids[id] = i + 1
      push!(data.D.V, Vertex(id, parse(Float64, parts[2]), parse(Float64, parts[3])))
    end
    # get blocks
    readline(f)
    for i in 1:nBlocks
      parts = split(readline(f), [',', ' ']; limit=0, keepempty=false)
      block = Vi([ids[parse(Int, part)] for part in parts[begin:end - 1]])
      push!(data.B, block)
      # define profit
      data.profits[block] = parse(Float64, parts[end])
    end
  end
  # update arcs ids
  data.D.distance = Dict{Tuple{Int, Int}, Float64}((ids[a[1]], ids[a[2]]) => distance for (a, distance) in data.D.distance)
  data.D.A = collect(keys(data.D.distance))
  # compact
  Vb = Set{Int}([i for b in data.B for i in b])
  data′, paths = compact(data, Vb)
  # update arcs
  data.D.A = collect(Set{Tuple{Int, Int}}(vcat(
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
  # return
  return data, Dict{Int, Int}(v => k for (k, v) in ids), data′, paths
end

function check_sbrp_complete_feasibility(data_complete::SBRPData, Vb)
  A′ = keys(data_complete.D.distance)
  return all((i, j) in A′ && (j, i) in A′ for i in Vb for j in Vb if i < j)
end

end
