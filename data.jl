import Unicode

mutable struct Vertex
   id_vertex::Int64
   pos_x::Float64
   pos_y::Float64
end

# Directed graph
mutable struct InputDigraph
   V::Array{Vertex} # set of vertices
   A::Array{Tuple{Int64,Int64}} # set of arcs
   distance::Dict{Tuple{Int64, Int64}, Float64} # time in 40 KM/h
end

δ⁺(A::Array{Tuple{Int64, Int64}}, i::Int64) = [(j, k) for (j, k) in A if j == i]
δ⁺(A::Array{Tuple{Int64, Int64}}, S::Set{Int64}) = [(i, j) for (i, j) in A if i in S && !in(j, S)]
δ⁺(A::Array{Tuple{Int64, Int64}}, S::Array{Int64}) = [(i, j) for (i, j) in A if i in S && !in(j, S)]
δ⁻(A::Array{Tuple{Int64, Int64}}, i::Int64) = [(j, k) for (j, k) in A if k == i]
δ⁻(A::Array{Tuple{Int64, Int64}}, S::Set{Int64}) = [(j, k) for (j, k) in A if k == i]
δ⁻(A::Array{Tuple{Int64, Int64}}, S::Array{Int64}) = [(j, k) for (j, k) in A if k == i]
time(data, a) = data.D.distance[a] / 40.0
time_block(data, block) = 4 * (sum(time(data, (block[i - 1], block[i])) for i in 2:length(block)) + time(data, (block[end], block[1])))

mutable struct SBRPData
  D::InputDigraph
  depot::Int64
  B::Array{Array{Int64, 1}}
  T::Float64
end

function compact(data::SBRPData)
  Vb = [i for b in data.B for i in b]
  data′, V, paths = SBRPData(
    InputDigraph(
                 Array{Vertex, 1}(vcat(data.D.V[data.depot], [data.D.V[i] for i in Vb])), 
                 Array{Tuple{Int64,Int64}, 1}(), 
                 Dict{Tuple{Int64, Int64}, Float64}()
                ), 
    data.depot, data.B, data.T # 2 hours
  ), 1:length(data.D.V), Dict{Tuple{Int64, Int64}, Array{Int64}}()
  for i in Vb
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
    for j in Vb
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
  [(data′.D.distance[(data.depot, i)] = data′.D.distance[(i, data.depot)] = 0.0) for i in Vb]
  data′.D.A = collect(keys(data′.D.distance))
  return data′, paths
end

function readSBRPData(app::Dict{String,Any}, compactPaths = false)
  depot = 1
  data, blocks, ids = SBRPData(
    InputDigraph(
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
      for j in 5:length(parts)
        id_block = parse(Int64, parts[j])
        !haskey(blocks, id_block) && (blocks[id_block] = Array{Tuple{Int64, Int64}, 1}())
        push!(blocks[id_block], a)
      end
    end
  end
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
  end
  # dummy arcs
  [data.D.distance[(depot, i)] = data.D.distance[(i, data.depot)] = 0.0 for i in 2:length(data.D.V)]
  # add arcs
  data.D.A = [keys(data.D.distance)...]
  # ids
  ids′ = Dict{Int64, Int64}(v => k for (k, v) in ids)
  # compact
  compactPaths && return data, ids′, compact(data)...
  # return
  return data, ids′
end

