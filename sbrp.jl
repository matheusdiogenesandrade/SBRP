module SBRP

include("symbols.jl")

using ..Data
#import ..InputDigraph

using BrkgaMpIpr

export time_block, SBRPData, compact, readSBRPData, readSBRPDataCarlos, distance_block, tour_distance, tour_time

mutable struct SBRPData <: AbstractInstance
  D::Data.InputDigraph
  depot::Int
  B::Vector{Vi}
  T::Float64
  profits::Dict{Vi, Float64}
end

# 40 km / 60 min = 40 km/h
NORMAL_SPEED = (40.0 * 1e3)/60.0

# arc time in 40 km/h
time(data, a) = a[1] == a[2] ? 0.0 : data.D.distance[a] / NORMAL_SPEED

# block distance in meters
distance_block(data::SBRPData, block::Vi) = sum(data.D.distance[(block[i - 1], block[i])] for i in 2:length(block)) + data.D.distance[(block[end], block[begin])]

# block service time made in 10 km/h
time_block(data::SBRPData, block::Vi) = 4 * distance_block(data, block) / NORMAL_SPEED 

# distance of a tour
tour_distance(data::SBRPData, tour::Vi) = sum(data.D.distance[(tour[i - 1], tour[i])] for i in 2:length(tour))

# tour time considering the time to spraying all the blocks
tour_time(data::SBRPData, tour::Vi, B::Vector{Vi}) = ((tour_distance(data, tour) - sum(distance_block(data, block) for block in B)) / NORMAL_SPEED) + sum(time_block(data, block) for block in B)

# get nodes belonging to some block
get_blocks_nodes(data::SBRPData) = Si([i for block in data.B for i in block])

# add dummy arcs connected to the nodes belonging to some block
add_dummy_arcs(data::SBRPData) = [data.D.distance[(data.depot, i)] = data.D.distance[(i, data.depot)] = 0.0 for i in get_blocks_nodes(data)]

function create_dummy_digraph(data″::SBRPData, paths″::Dict{Arc, Vi})

    V = collect(data″.D.V)

    new_id = max()

end

function create_no_one_degree_paths_digraph(data::SBRPData)

  # setup
  V, A, Vb, distance = data.D.V, ArcsSet(data.D.A), get_blocks_nodes(data), data.D.distance

  # utils
  valid_node(A::ArcsSet, i::Int) = ∧(∨(∧(v⁻(A, i) ⊆ v⁺(A, i), d⁺(A, i) <= 2, d⁻(A, i) <= 2), ∧(v⁺(A, i) ⊆ v⁻(A, i), d⁺(A, i) <= 2, d⁻(A, i) <= 2)), !in(i, Vb))

  # new digraph by selecting the proper nodes
  V′ = Si([i for i in keys(V) if valid_node(A, i)])
  A′ = filter(a -> a[1] in V′ && a[2] in V′, A)
  
  # get paths
  paths, visited_nodes = Dict{Arc, Vi}(), Si()
  function dfs(curr::Int, path::Vi)
    # base case
    !valid_node(A′, curr) && error("Node $curr is not valid")

    # children
    for next in v⁺(A′, curr)
      next in visited_nodes && continue 

      # store node
      push!(path, next)
      push!(visited_nodes, next)

      # go on
      dfs(next, path)
    end
  end
  for i in V′
    # checking if it is not a head node
    !(d⁺(A′, i) == 1 && ⊆(v⁻(A′, i), v⁺(A′, i))) && continue

    # dfs
    path = Vi([i])
    push!(visited_nodes, i)
    dfs(i, path)

    # store path
    if length(path) >= 2 
      paths[(path[begin], path[end])] = path

      # check if reversed is feasible
      path′ = reverse(path)
      all((path′[i], path′[i + 1]) in A for i in 1:length(path′) - 1) && (paths[(path′[begin], path′[end])] = path′)
    end
  end

  # connect with Vb's nodes
  function in_node(a)
    arr = setdiff(v⁻(A, a[1]), v⁻(A′, a[1]))
    return length(arr) == 1 ? first(arr) : nothing
  end
  function out_node(a)
    arr = setdiff(v⁺(A, a[2]), v⁺(A′, a[2]))
    return length(arr) == 1 ? first(arr) : nothing
  end
  paths′ = filter(pair -> in_node(pair[1]) == out_node(pair[1]) == nothing, paths)
  paths″ = Dict{Arc, Vi}()
  for (a, path) in filter(pair -> !in(pair[1], keys(paths′)), paths)
    a′, path′, nodeIn, nodeOut = a, copy(path), in_node(a), out_node(a)
    if nodeIn != nothing
      a′ = (nodeIn, a′[2])
      path′ = ∪(nodeIn, path′)
    end
    if nodeOut != nothing
      a′ = (a′[1], nodeOut)
      path′ = ∪(path′, nodeOut)
    end
    paths″[a′] = path′
  end
  paths = merge(paths′, paths″)

  # build new digraph
  data″ = deepcopy(data)
  A″ = collect(∪(filter(a -> !(a[1] in visited_nodes && a[2] in visited_nodes), A), keys(paths)))
#  A″ = collect(A)
  D″ = Data.InputDigraph(
    Dict{Int, Vertex}(i => V[i] for i in vcat([i for (i, j) in A″], [j for (i, j) in A″])),
    A″,
    ArcCostMap(a => a in keys(paths) ? ∑(distance[(paths[a][i - 1], paths[a][i])] for i in 2:length(paths[a])) : distance[a] for a in A″)
  )
  data″.D = D″

  # remove begins and endings
  [paths[a] = path[begin + 1 : end - 1] for (a, path) in paths]

  return data″, paths
end

function create_complete_digraph(data::SBRPData)

  # setup
  A, B = data.D.A, data.B
  Vb, depot = Si([i for b in B for i in b]), data.depot

  data′, V, paths = SBRPData(
    Data.InputDigraph(
                 Dict{Int, Vertex}(vcat((depot => data.D.V[depot]), [(i => data.D.V[i]) for i in Vb])...), 
                 Arcs(), 
                 ArcCostMap((i, j) => 0.0 for (i, j) in χ(Vb) if i != j)
                ), 
    depot, 
    B, 
    data.T, # 2 hours
    data.profits
  ), keys(data.D.V), Dict{Arc, Vi}((i, j) => [] for (i, j) in χ(Vb) if i != j)

  # adjList
  adjList = Dict{Int, Vi}(i => Vi() for i in keys(data.D.V))
  for (i, j) in A
    push!(adjList[i], j)
  end

  adjList = Dict{Int64, Vi}(i => [] for i in V)
  [push!(adjList[i], j) for (i, j) in A]

  # get paths
  for i in Vb
#	flush_println(i)
    # bfs
    distances, pred, q = Dict{Int, Float64}(i => typemax(Float64) for i in V), Dict{Int, Int}(i => i for i in V), [i]
    distances[i] = 0.0
    while !∅(q)
      curr = popfirst!(q)
  #    for (curr, next) in δ⁺(A, curr)
      for next in adjList[curr]
#        if !in(next, [data.depot, i]) && (pred[next] == next || distances[next] > distances[curr] + data.D.distance[(curr, next)])
        if !in(next, [data.depot, i]) && distances[next] > distances[curr] + data.D.distance[(curr, next)]
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

  # dummy arcs
  [data′.D.distance[(i, depot)] = data′.D.distance[(depot, i)] = 0.0 for i in Vb]

  # add arcs
  data′.D.A = collect(keys(data′.D.distance))

  # remove beginnings and endings of each path
  for (a, path) in paths
    popfirst!(path)
    pop!(path)
  end

  # check feasibility
  check_sbrp_complete_feasibility(data′, Vb)

  return data′, paths
end

function readSBRPDataCarlos(app::Dict{String,Any})
  data, blocks = SBRPData(
    Data.InputDigraph(
                 Dict{Int, Vertex}(1 => Vertex(-1, -1, -1)), 
                 Arcs(), 
                 ArcCostMap()
                ), 
    -1,
    VVi(),
    parse(Float64, app["vehicle-time-limit"]),
    ArcCostMap()
  ), Dict{Int, Arcs}()

  open(app["instance"]) do f

    nNodes::Int, nArcs::Int, nBlocks::Int = map(str -> parse(Int, str), split(readline(f), [' ']; limit=3, keepempty=false))

    flush_println("Reading nodes")
    # get vertices
    for i in 1:nNodes

        strs = split(readline(f), [' ']; keepempty=false)

        id = parse(Int, strs[2])
        x = parse(Float64, strs[3])
        y = parse(Float64, strs[4])

        data.D.V[id] = Vertex(id, x, y)
    end

    flush_println("Reading arcs")
    # get arcs
    block_arcs = Dict{Int, Vector{Tuple{Int, Int}}}()
    block_profit = Dict{Int, Int}()

    for n in 1:nArcs

        strs = split(readline(f), [' ']; keepempty=false)

        a = (parse(Int, strs[2]), parse(Int, strs[3]))
        distance = parse(Float64, strs[4])
        block_id = parse(Int, strs[5])
        cases = parse(Int, strs[6])

        data.D.distance[a] = distance

        block_id == -1 && continue

        # new block
        if !haskey(block_profit, block_id)
            block_profit[block_id] = 0
            block_arcs[block_id] = Vector{Tuple{Int, Int}}()
        end

        block_profit[block_id] += cases

        push!(block_arcs[block_id], a)

    end

    # get blocks
    function get_block(arcs::Vector{Tuple{Int, Int}})
        block = Vi()
        curr, next = first(arcs)

        push!(block, curr)
        deleteat!(arcs, 1)

        while true
            idx = findfirst((i, j)::Tuple{Int, Int} -> next == i, arcs)
            idx == nothing && break

            curr, next = arcs[idx]

            push!(block, curr)
            deleteat!(arcs, idx)
        end

        return block
    end

    for (block_id, arcs) in block_arcs

        block = get_block(arcs)
        push!(data.B, block)

        data.profits[block] = block_profit[block_id]
    end

  end

  # update arcs ids
  data.D.A = collect(keys(data.D.distance))

  # add dummy depot
  depot = data.depot = maximum(collect(keys(data.D.V))) + 1
  data.D.V[depot] = Vertex(depot, -1, -1)

  # dummy weights
  add_dummy_arcs(data)

  flush_println("Checking feasibility")
  # check feasibility
  Vb = get_blocks_nodes(data)
  distances = calculate_shortest_paths(data)
  any(!in((i, j), keys(distances)) for (i, j) in χ(Vb)) && error("The SBRP instance it is not connected")

  flush_println("Checking complete graph")
  # compact in complete graph
  data′, paths′ = create_complete_digraph(data)

  # check feasibility
  check_feasibility(data, data′)

  # consider only shortest paths arcs
    # copy
#  data‴ = deepcopy(data)
    # set arcs
#  data‴.D.A = collect(ArcsSet(vcat(
#                                   [(path[i], path[i + 1]) for (a, path) in paths′ for i in 1:(length(path) - 1)], # min paths arcs
#                                   [(a[1], path[begin]) for (a, path) in paths′ if !∅(path)], # min paths arcs
#                                   [(path[end], a[2]) for (a, path) in paths′ if !∅(path)], # min paths arcs
#                                   [a for (a, path) in paths′ if ∅(path)], # min paths arcs (edge case)
#                                   [(b[i], b[i + 1]) for b in data.B for i in 1:(length(b) - 1)], # blocks arcs
#                                   [(b[end], b[begin]) for b in data.B], # blocks arcs
#                                   [(data.depot, i) for i in Vb], # depot arcs
#                                   [(i, data.depot) for i in Vb] # depot arcs
#                                  )))
#  data‴.D.distance = ArcCostMap(a => data.D.distance[a] for a in data‴.D.A if a[1] != a[2])

  # check feasibility
#  check_feasibility(data, data‴)

  # set
#  data = data‴

  # compact algorithm
#  data″, paths″ = create_no_one_degree_paths_digraph(data)

  # check feasibility
#  check_feasibility(data, data″)

  # return
#  return data, data′, paths′, data″, paths″
  return data, data′, paths′
end

function readSBRPData(app::Dict{String,Any})
  depot = 1
  data, blocks = SBRPData(
    Data.InputDigraph(
                 Dict{Int, Vertex}(1 => Vertex(-1, -1, -1)), 
                 Arcs(), 
                 ArcCostMap()
                ), 
    depot,
    VVi(),
    parse(Float64, app["vehicle-time-limit"]),
    ArcCostMap()
  ), Dict{Int, Arcs}()

  open(app["instance"]) do f

    # ignore the first two lines
    readline(f)
    readline(f)
    nNodes, nArcs, nBlocks = parse(Int, split(readline(f), [' ']; limit=0, keepempty=false)[end]), parse(Int, split(readline(f), [' ']; limit=0, keepempty=false)[end]) + parse(Int, split(readline(f), [' ']; limit=0, keepempty=false)[end]), parse(Int, split(readline(f), [' ']; limit=0, keepempty=false)[end])
    readline(f)
    readline(f)

    flush_println("Reading arcs")
    # get arcs
    for n in 1:nArcs
      parts = split(readline(f), ['(', ')', ' ', ',']; limit=0, keepempty=false)
      data.D.distance[(parse(Int, parts[1]), parse(Int, parts[2]))] = floor(Int, parse(Float64, parts[end]))
    end

    flush_println("Reading nodes")
    # get vertices
    readline(f)
    for i in 1:nNodes
      parts = split(readline(f), [' ']; limit=0, keepempty=false)
      id = parse(Int, parts[1])
      data.D.V[id] = Vertex(id, parse(Float64, parts[2]), parse(Float64, parts[3]))
    end

    flush_println("Reading blocks")
    # get blocks
    readline(f)
    for i in 1:nBlocks
      parts = split(readline(f), [',', ' ']; limit=0, keepempty=false)
      block = Vi([parse(Int, part) for part in parts[begin:end - 1]])
      push!(data.B, block)

      # define profit
      data.profits[block] = app["unitary-profits"] ? 1.0 : parse(Float64, parts[end])
      
    end

  end

  # update arcs ids
  data.D.A = collect(keys(data.D.distance))

  # dummy weights
  add_dummy_arcs(data)

  flush_println("Checking feasibility")
  # check feasibility
  Vb = get_blocks_nodes(data)
  distances = calculate_shortest_paths(data)
  any(!in((i, j), keys(distances)) for (i, j) in χ(Vb)) && error("The SBRP instance it is not connected")

  flush_println("Checking complete graph")
  # compact in complete graph
  data′, paths′ = create_complete_digraph(data)

  # check feasibility
  check_feasibility(data, data′)

  # consider only shortest paths arcs
    # copy
  data‴ = deepcopy(data)
    # set arcs
  data‴.D.A = collect(ArcsSet(vcat(
                                   [(path[i], path[i + 1]) for (a, path) in paths′ for i in 1:(length(path) - 1)], # min paths arcs
                                   [(a[1], path[begin]) for (a, path) in paths′ if !∅(path)], # min paths arcs
                                   [(path[end], a[2]) for (a, path) in paths′ if !∅(path)], # min paths arcs
                                   [a for (a, path) in paths′ if ∅(path)], # min paths arcs (edge case)
                                   [(b[i], b[i + 1]) for b in data.B for i in 1:(length(b) - 1)], # blocks arcs
                                   [(b[end], b[begin]) for b in data.B], # blocks arcs
                                   [(depot, i) for i in Vb], # depot arcs
                                   [(i, depot) for i in Vb] # depot arcs
                                  )))
  data‴.D.distance = ArcCostMap(a => data.D.distance[a] for a in data‴.D.A)

  # check feasibility
  check_feasibility(data, data‴)

  # set
  data = data‴

  # compact algorithm
#  data″, paths″ = create_no_one_degree_paths_digraph(data)

  # check feasibility
#  check_feasibility(data, data″)

  # return
#  return data, data′, paths′, data″, paths″
  return data, data′, paths′
end

function check_feasibility(data::SBRPData, data′::SBRPData) 
  Vb = get_blocks_nodes(data)
  distances = calculate_shortest_paths(data)
  distances′ = calculate_shortest_paths(data′)
  any(!in((i, j), keys(distances′)) for (i, j) in χ(Vb)) && error("The SBRP instance it is not connected")
  for a in χ(Vb)
      abs(distances[a] - distances′[a]) > 1e-3 && error("Original digraph distance is $(distances[a]) while the new digraph distance is $(distances′[a])")
  end
end

function calculate_shortest_paths(data::SBRPData)
  V, A, Vb = data.D.V, data.D.A, get_blocks_nodes(data)

  adjList = Dict{Int64, Vi}(i => [] for i in keys(V))
  [push!(adjList[i], j) for (i, j) in A]

  distances = Dict{Arc, Float64}()
  for i in Vb
#  	flush_println(i)
    # bfs
    q = [i]
    distances[(i, i)] = 0.0
    while !∅(q)
      # curr node
      curr = pop!(q)
      # next nodes
      for next in adjList[curr] 
        # edge case
        next == data.depot && continue

        # update distance and push new nodes
        b = distances[(i, curr)] + data.D.distance[(curr, next)]
        if !in((i, next), keys(distances)) || distances[(i, next)] > b
          distances[(i, next)] = b
          push!(q, next)
        end
      end
    end
  end

  return distances
end

check_sbrp_complete_feasibility(data::SBRPData, Vb) = any([!in((i, j), data.D.A) for i in Vb for j in Vb if i != j]) ? error("The Complete SBRP instance is not feasible") : ()

function check_sbrp_feasibility(data::SBRPData)
  V′ = get_blocks_nodes(data)
  A = data.D.A
  # Kosaraju’s algorithm
  # outgoing arcs
  # bfs
  v = first(V′)
  mandatory_visited_vertices = Si([v])
  visited_vertices = Si([v])
  q = [v]
  while !∅(q)
    curr = popfirst!(q) 
    for (i, j) in δ⁺(A, curr)
      if !in(j, visited_vertices)
        push!(q, j)
        push!(visited_vertices, j)
        if j in V′
          push!(mandatory_visited_vertices, j)
          if length(mandatory_visited_vertices) == length(V′)
            empty(q)
            break
          end
        end
      end
    end
  end
  if length(mandatory_visited_vertices) < length(V′)
    for u in V′
      !in(u, mandatory_visited_vertices) && error("There is not path between $v and $u")
    end
  end
  # Kosaraju’s algorithm
  # ingoing arcs
  # bfs
  mandatory_visited_vertices = Si([v])
  visited_vertices = Si([v])
  q = [v]
  while !∅(q)
    curr = popfirst!(q) 
    for (i, j) in δ⁻(A, curr)
      if !in(i, visited_vertices)
        push!(q, i)
        push!(visited_vertices, i)
        if i in V′
          push!(mandatory_visited_vertices, i)
          if length(mandatory_visited_vertices) == length(V′)
            empty(q)
            break
          end
        end
      end
    end
  end
  if length(mandatory_visited_vertices) < length(V′)
    for u in V′
      !in(u, mandatory_visited_vertices) && error("There is not path between $v and $u")
    end
  end
end

end
