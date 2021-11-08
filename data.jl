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
δ⁺(A::Array{Tuple{Int64, Int64}}, S::Set{Int64}) = [(j, k) for (j, k) in A if j == i]
δ⁺(A::Array{Tuple{Int64, Int64}}, S::Array{Int64}) = [(j, k) for (j, k) in A if j == i]
δ⁻(A::Array{Tuple{Int64, Int64}}, i:Int64) = [(j, k) for (j, k) in A if k == i]
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

function readSBRPData(app::Dict{String,Any})
  depot = 1
  data = SBRPData(
    InputDigraph(
                 Array{Vertex, 1}([Vertex(-1, -1, -1)]), 
                 Array{Tuple{Int64,Int64}, 1}(), 
                 Dict{Tuple{Int64, Int64}, Float64}()
                ), 
    -1,
    Array{Array{Int64, 1}, 1}(),
    120.0 # 2 hours
  )
  blocks = Dict{Int64, Array{Tuple{Int64, Int64}, 1}}()
#  Vₘ = Dict{Int64, Int64}()
  open(app["instance"]) do f
    parts = split(readline(f), [' ']; limit=0, keepempty=false)
    # get params
    nNodes, nArcs, nBlocks = parse(Int64, parts[1]), parse(Int64, parts[2]), parse(Int64, parts[3])
    # get nodes
    data.D.V = vcat(data.D.V, [Vertex(parse(Int64, split(readline(f), [' ']; limit=0, keepempty=false)[2]), 0.0, 0.0) for i in 1:nNodes])
    # get distances
    for i in 1:nArcs
      parts = split(readline(f), [' ']; limit=0, keepempty=false)
      a = (parse(Int64, parts[2]), parse(Int64, parts[3]))
      data.D.distance[a] = parse(Float64, parts[4])
      # get blocks
      if parts[5] != "-1"
        for j in 5:length(parts)
          id_block = parse(Int64, parts[j])
          if !haskey(blocks, id_block)
            blocks[id_block] = Array{Tuple{Int64, Int64}, 1}()
          end
          push!(blocks[id_block], a)
        end
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
  # dommy arcs
  for v in data.D.V
    data.D.distance[(data.depot, v.id_vertex)] = data.D.distance[(v.id_vertex, data.depot)] = 0.0
  end
  # add arcs
  data.D.A = [keys(A)...]
  # return
  return data
end

