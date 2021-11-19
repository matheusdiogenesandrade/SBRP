module SBRP

using ..Data
#import ..InputDigraph

export time_block, SBRPData, compact, readSBRPData, checkSBRPfeasibility

time_block(data, block) = 4 * (sum(Data.time(data, (block[i - 1], block[i])) for i in 2:length(block)) + Data.time(data, (block[end], block[1])))

mutable struct SBRPData
  D::Data.InputDigraph
  depot::Int64
  B::Array{Array{Int64, 1}}
  T::Float64
end

function compact(data::SBRPData, V′)
  data′, V, paths = SBRPData(
    Data.InputDigraph(
                 Array{Vertex, 1}(vcat(data.D.V[data.depot], [data.D.V[i] for i in V′])), 
                 Array{Tuple{Int64,Int64}, 1}(), 
                 Dict{Tuple{Int64, Int64}, Float64}()
                ), 
    data.depot, data.B, data.T # 2 hours
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
  return data′, paths
end

function readSBRPData(app::Dict{String,Any})
  depot = 1
  data, blocks, ids = SBRPData(
    Data.InputDigraph(
                 Array{Vertex, 1}([Vertex(-1, -1, -1)]), 
                 Array{Tuple{Int64,Int64}, 1}(), 
                 Dict{Tuple{Int64, Int64}, Float64}()
                ), 
    depot,
    Array{Array{Int64, 1}, 1}(),
    120.0 # 2 hours
  ), Dict{Int64, Array{Tuple{Int64, Int64}, 1}}(), Dict{Int64, Int64}()
#  Vₘ = Dict{Int64, Int64}()
  open(app["instance"]) do f
    parts = split(readline(f), [' ']; limit=0, keepempty=false)
    # get params
    nNodes, nArcs, nBlocks = parse(Int64, parts[1]), parse(Int64, parts[2]), parse(Int64, parts[3])
    # get nodes
    for i in 1:nNodes
      id = parse(Int64, split(readline(f), [' ']; limit=0, keepempty=false)[2])
      ids[id] = i + 1
      push!(data.D.V, Vertex(id, 0.0, 0.0))
    end
    # get distances
    for i in 1:nArcs
      parts = split(readline(f), [' ']; limit=0, keepempty=false)
      a = (ids[parse(Int64, parts[2])], ids[parse(Int64, parts[3])])
      data.D.distance[a] = parse(Float64, parts[4])
      # get blocks
      parts[5] == "-1" && continue
      id_block = parse(Int64, parts[5])
      !haskey(blocks, id_block) && (blocks[id_block] = Array{Tuple{Int64, Int64}, 1}())
      push!(blocks[id_block], a)
    end
  end
  n_blocks, k = 24, 1
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

    k >= n_blocks && break
    k = k + 1
  end
  # add arcs
  data.D.A = [keys(data.D.distance)...]
  # compact
  Vb = Set{Int64}([i for b in data.B for i in b])
  data′, paths = compact(data, Vb)
  # update arcs
  data.D.A = vcat(
    collect(Set{Tuple{Int64, Int64}}([(path[i], path[i + 1]) for (a, path) in paths for i in 1:(length(path) - 1)])),
    [(depot, i) for i in Vb],
    [(i, depot) for i in Vb]
  )
  # dummy weights
  [data.D.distance[(depot, i)] = data.D.distance[(i, depot)] = 0.0 for i in Vb]
  # check feasibility
  !checkSBRPfeasibility(data′, Vb) && error("The SBRP instance is not feasible")
  # return
  return data, Dict{Int64, Int64}(v => k for (k, v) in ids), data′, paths
end

function checkSBRPfeasibility(data_complete::SBRPData, V′)
  A′ = keys(data_complete.D.distance)
  return all((i, j) in A′ && (j, i) in A′ for i in V′ for j in V′ if i < j)
end

end
