using BrkgaMpIpr

mutable struct SBRPData <: AbstractInstance
    D::InputDigraph
    depot::Int
    B::Vector{Vi}
    T::Float64
    profits::Dict{Vi, Float64}
end

mutable struct SBRPSolution # SBRP Solution
   tour::Vi                 # route
   B::VVi                   # serviced blocks
end

# 40 km / 60 min = 40 km/h
const NORMAL_SPEED::Float64 = (40.0 * 1e3)/60.0

# arc time in 40 km/h
const time(data::SBRPData, (i, j)::Arc)::Float64 = i == j ? 0.0 : data.D.distance[Arc(i, j)] / NORMAL_SPEED

# arc distance
const distance(data::SBRPData, (i, j)::Arc)::Float64 = i == j ? 0.0 : data.D.distance[Arc(i, j)]

# arc time with custom speed
#const time(data::SBRPData, (i, j)::Arc, speed::Float64)::Float64 = i == j ? 0.0 : data.D.distance[Arc(i, j)] / speed

# block distance in meters
const blockDistance(data::SBRPData, block::Vi)::Float64 = sum(a::Arc -> data.D.distance[a], Arcs(collect(zip(block[begin:end - 1], block[begin + 1:end])))) + data.D.distance[Arc(last(block), first(block))]

# block service time made in 10 km/h
const blockTime(data::SBRPData, block::Vi)::Float64 = 4 * blockDistance(data, block) / NORMAL_SPEED 

# distance of a tour
const tourDistance(data::SBRPData, tour::Vi)::Float64 = sum(a::Arc -> data.D.distance[a], Arcs(collect(zip(tour[begin:end - 1], tour[begin + 1:end]))))

# tour time considering the time to spraying all the blocks
const tourTime(data::SBRPData, solution::SBRPSolution)::Float64 = ((tourDistance(data, solution.tour) - sum(block::Vi -> blockDistance(data, block), solution.B)) / NORMAL_SPEED) + sum(block::Vi -> blockTime(data, block), solution.B)

# get nodes belonging to some block
const getBlocksNodes(data::SBRPData)::Si = Si(reduce(vcat, data.B))

# add dummy arcs connected to the nodes belonging to some block
function addDummyArcs(data::SBRPData) 

    @debug "Adding dummy arcs"

    for i::Int in getBlocksNodes(data)

        @debug ">> For node $i"

        data.D.distance[Arc(data.depot, i)] = data.D.distance[Arc(i, data.depot)] = 0.0 
    end
end

function createDummyDigraph(data″::SBRPData, paths″::Dict{Arc, Vi})

    @debug "Creating dummy digraph"

    V::Vi = collect(data″.D.V)

    #    new_id = max()

end

#=
Create a SBRP instance in which there is no path in which all the nodes have degree (in or out) equals to one
input:
- data::SBRPData is a SBRP instance
output:
- data″::SBRPData is the new SBRP instance
- paths::Dict{Pair{Int, Int}, Vector{Int}} is the relation of paths
=#
function createNoOneDegreePathsDigraph(data::SBRPData)

    @debug "Creating no one degree paths digraph"

    # setup
    V::Dict{Int, Vertex} = data.D.V
    A::ArcsSet = ArcsSet(data.D.A)
    Vb::Si = getBlocksNodes(data)
    distance::ArcCostMap = data.D.distance

    # utils
    v⁺::Dict{Int, Si} = Dict{Int, Si}(map(i::Int -> i => Si(), V))
    v⁻::Dict{Int, Si} = Dict{Int, Si}(map(j::Int -> j => Si(), V))
    d⁺(i::Int)::Int = length(v⁺[i])
    d⁻(j::Int)::Int = length(v⁻[j])

    for (i::Int, j::Int) in A
        push!(v⁺[i], j)
        push!(v⁻[j], i)
    end

    validNode(i::Int)::Bool = ∧(
                                ∨(
                                  ∧(
                                    v⁻[i] ⊆ v⁺[i], 
                                    d⁺(i) <= 2, 
                                    d⁻(i) <= 2
                                   ), 
                                  ∧(
                                    v⁺[i] ⊆ v⁻[i], 
                                    d⁺(i) <= 2, 
                                    d⁻(i) <= 2
                                   )
                                 ), 
                                !in(i, Vb)
                               )

    # new digraph by selecting the proper nodes
    V′::Si = Si(filter(i::Int -> validNode(A, i), keys(V)))
    A′::ArcsSet = ArcsSet(filter((i, j)::Arc -> i in V′ && j in V′, A))

    v⁺′::Dict{Int, Si} = Dict{Int, Si}(map(i::Int -> i => Si(), V′))
    v⁻′::Dict{Int, Si} = Dict{Int, Si}(map(j::Int -> j => Si(), V′))
    d⁺′(i::Int)::Int = length(v⁺′[i])
    d⁻′(j::Int)::Int = length(v⁻′[j])

    for (i::Int, j:Int) in A′
        push!(v⁺′[i], j)
        push!(v⁻′[j], i)
    end

    validNode′(i::Int)::Bool = ∧(
                                ∨(
                                  ∧(
                                    v⁻′[i] ⊆ v⁺′[i], 
                                    d⁺′(i) <= 2, 
                                    d⁻′(i) <= 2
                                   ), 
                                  ∧(
                                    v⁺′[i] ⊆ v⁻′[i], 
                                    d⁺′(i) <= 2, 
                                    d⁻′(i) <= 2
                                   )
                                 ), 
                                !in(i, Vb)
                               )

    # get paths
    paths::Dict{Arc, Vi} = Dict{Arc, Vi}()
    visited_nodes::Si = Si()

    function dfs(curr::Int, path::Vi)

        # base case
        if !validNode′(curr)
            throw(InvalidStateException("Node $curr is not valid"))
        end

        # children
        for next::Int in v⁺′[curr]

            next in visited_nodes && continue 

            # store node
            push!(path, next)
            push!(visited_nodes, next)

            # go on
            dfs(next, path)

        end

    end

    for i::Int in V′

        @debug "Obtaining one degree path starting at node $i"

        # checking if it is not a head node
        !(d⁺′(i) == 1 && ⊆(v⁻′[i], v⁺′[i])) && continue

        # dfs
        path::Vi = Vi([i])
        push!(visited_nodes, i)
        dfs(i, path)

        # store path
        if length(path) >= 2 

            @debug ">> Saving path"

            paths[(path[begin], path[end])] = path

            # check if reversed is feasible
            path′::Vi = reverse(path)

            if all((i, j)::Arc -> (i, j) in A, zip(path′[begin:end - 1], path′[begin + 1:end]))

                paths[(path′[begin], path′[end])] = path′

            end

        end
    end

    # connect with Vb's nodes
    function inNode((i, j)::Arc)::Union{Int, Nothing}

        arr::Vi = setdiff(v⁻[i], v⁻′[i])

        return length(arr) == 1 ? first(arr) : nothing

    end

    function outNode((i, j)::Arc)::Union{Int, Nothing}

        arr::Vi = setdiff(v⁺[j], v⁺′[j])

        return length(arr) == 1 ? first(arr) : nothing

    end

    paths′::Dict{Arc, Vi} = filter((a, path)::Pair{Arc, Vi} -> inNode(a) == outNode(a) == nothing, paths)
    paths″::Dict{Arc, Vi} = Dict{Arc, Vi}()

    @debug "Mapping paths to arcs"
    for (a::Arc, path::Vi) in filter((a, _)::Pair{Arc, Vi} -> !in(a, keys(paths′)), paths)

        @debug ">> Mapping arc $a"

        a′::Arc = a
        path′::Vi = copy(path)
        nodeIn::Union{Int, Nothing} = inNode(a)
        nodeOut::Union{Int, Nothing} = outNode(a)

        if nodeIn != nothing
            a′ = (nodeIn, a′[2])
            pushfirst!(path′, nodeIn)
        end

        if nodeOut != nothing
            a′ = (a′[1], nodeOut)
            push!(path′, nodeOut)
        end

        paths″[a′] = path′

    end

    paths = merge(paths′, paths″)

    # build new digraph
    data″::SBRPData = deepcopy(data)
    A″::Arcs = collect(∪(filter((i, j)::Arc -> !(i in visited_nodes && j in visited_nodes), A), keys(paths)))

    D″ = InputDigraph(
                      Dict{Int, Vertex}(map(i::Int -> i => V[i], vcat(map((i, j)::Arc -> i, A″), map((i, j)::Arc -> j, A″)))),
                      A″,
                      ArcCostMap(map(a::Arc -> a => a in keys(paths) ? sum((i, j)::Arc -> distance[Arc(i, j)], zip(paths[a][begin:end - 1], paths[a][begin + 1:end])) : distance[a], A″))
                     )

    data″.D = D″

    # remove begins and endings
    for (a::Arc, path::Vi) in paths
        paths[a] = path[begin + 1 : end - 1]
    end

    return data″, paths
end

#=
Create a SBRP instance in which the digraph is complete and metric
input:
- data::SBRPData is a SBRP instance
output:
- data´::SBRPData is the new SBRP instance
- paths´::Dict{Pair{Int, Int}, Vector{Int}} is the relation of paths
=#
function createCompleteDigraph(data::SBRPData)::Pair{SBRPData, Dict{Arc, Vi}}

    @debug "Creating complete digraph"

    # setup
    A::Arcs = data.D.A
    arcs_set::ArcsSet = ArcsSet(data.D.A)
    B::VVi = data.B
    Vb::Si = Si([i for b in B for i in b])
    depot::Int = data.depot

    data′::SBRPData = SBRPData(
                               InputDigraph(
                                            Dict{Int, Vertex}(vcat(
                                                                   depot => data.D.V[depot], 
                                                                   map(i::Int -> i => data.D.V[i], collect(Vb))
                                                                  )), 
                                            Arcs(), 
                                            ArcCostMap(map(
                                                           a::Arc -> a => 0.0, 
                                                           filter((i, j)::Arc -> i != j, χ(Vb))
                                                          ))
                                           ), 
                               depot, 
                               B, 
                               data.T, # 2 hours
                               data.profits
                              )
    V::Vi = Vi(collect(keys(data.D.V)))
    paths::Dict{Arc, Vi} = Dict{Arc, Vi}(map((i, j)::Arc -> Arc(i, j) => Vi(), filter((i, j)::Pair{Int, Int} -> i != j, χ(Vb))))

    # adjList
    @debug "Creating adjacency list"

    adjList::Dict{Int, Vi} = Dict{Int, Vi}(map(i::Int -> i => Vi(), V))

    for (i::Int, j::Int) in A
        push!(adjList[i], j)
    end

    # subrotine to find all shortest paths with root at node i
    function findPaths(i::Int)

        @debug ">> Finding paths from node $i"

        # bfs
        distances::Dict{Int, Float64} = Dict{Int, Float64}(map(i::Int -> i => typemax(Float64), V))
        pred::Dict{Int, Int} = Dict{Int, Int}(map(i::Int -> i => i, V))
        q::Vi = Vi([i])

        distances[i] = 0.0

        while !isempty(q)

            curr::Int = popfirst!(q)

            for next::Int in adjList[curr]

                if !in(next, [data.depot, i]) && distances[next] > distances[curr] + data.D.distance[Arc(curr, next)]

                    distances[next] = distances[curr] + data.D.distance[Arc(curr, next)]
                    pred[next] = curr
                    push!(q, next)

                end

            end

        end

        # add cost and distance
        for j::Int in Vb

            pred[j] == j && continue

            data′.D.distance[Arc(i, j)] = distances[j]
            paths[Arc(i, j)] = Vi()
            curr::Int = j

            while pred[curr] != curr

                pushfirst!(paths[Arc(i, j)], curr)
                curr = pred[curr]

            end

            pushfirst!(paths[Arc(i, j)], i)

        end
    end

    # get paths in parallel
    @debug "Finding paths in parallel"

    path_tasks::Vector{Task} = Vector{Task}()
    for i::Int in Vb
        t = @task findPaths(i)
        schedule(t)
        push!(path_tasks, t)
    end

    # run tasks
    for t::Task in path_tasks
        wait(t)
    end

    # check feasibility of the generated paths
    for (_, path::Vi) in paths
        
        path_arcs::Vector{Tuple{Int, Int}} = collect(zip(path[begin:end - 1], path[begin + 1:end]))

        idx::Union{Nothing, Int} = findfirst((i, j)::Tuple{Int, Int} -> !in(Arc(i, j), arcs_set), path_arcs)

        if idx != nothing
            throw(InvalidStateException("In path $path, the arc $((path[idx], path[idx + 1])) does not exist", :closed))
        end

    end

    # dummy arcs
    @debug "Adding dummy arc weights"

    for i::Int in Vb
        data′.D.distance[Arc(i, depot)] = data′.D.distance[Arc(depot, i)] = 0.0
    end

    # add arcs
    data′.D.A = collect(keys(data′.D.distance))

    # remove beginnings and endings of each path
    for (a::Arc, path::Vi) in paths

        popfirst!(path)
        pop!(path)

    end

    # check feasibility
    checkSBRPCompleteFeasibility(data′, Vb)


    return Pair{SBRPData, Dict{Arc, Vi}}(data′, paths)
end

#=
Create a SBRP instance from a Carlo's instance
input:
- app::Dict{String, Any} is the paramenters relation
output:
- data::SBRPData is the SBRP instance
=#
function readSBRPDataCarlos(app::Dict{String, Any})::SBRPData

    @debug "Reading Carlos SBRP instances"

    data::SBRPData = SBRPData(
                              InputDigraph(
                                           Dict{Int, Vertex}(1 => Vertex(-1, -1, -1)), 
                                           Arcs(), 
                                           ArcCostMap()
                                          ), 
                              -1,
                              VVi(),
                              parse(Float64, app["vehicle-time-limit"]),
                              ArcCostMap()
                             )
    blocks::Dict{Int, Arcs} = Dict{Int, Arcs}()

    open(app["instance"]) do f::IOStream

        nNodes::Int, nArcs::Int, nBlocks::Int = map(str::String -> parse(Int, str), Vector{String}(split(readline(f), [' ']; limit=3, keepempty=false)))

        # get nodes
        @debug "Reading nodes"

        for i::Int in 1:nNodes

            strs::Vector{String} = split(readline(f), [' ']; keepempty=false)

            id::Int = parse(Int, strs[2])
            x::Float64 = parse(Float64, strs[3])
            y::Float64 = parse(Float64, strs[4])

            data.D.V[id] = Vertex(id, x, y)
        end

        # get arcs
        @debug "Reading arcs"

        block_arcs::Dict{Int, Arcs} = Dict{Int, Arcs}()
        block_profit::Dict{Int, Int} = Dict{Int, Int}()

        for n::Int in 1:nArcs

            strs::Vector{String} = split(readline(f), [' ']; keepempty=false)

            a::Arc = Arc(parse(Int, strs[2]), parse(Int, strs[3]))
            distance::Float64 = parse(Float64, strs[4])
            block_id::Int = parse(Int, strs[5])
            cases::Int = parse(Int, strs[6])

            data.D.distance[a] = distance

            block_id == -1 && continue

            # new block
            if !haskey(block_profit, block_id)
                block_profit[block_id] = 0
                block_arcs[block_id] = Arcs()
            end

            block_profit[block_id] += cases

            push!(block_arcs[block_id], a)

        end

        # get blocks
        @debug "Reading blocks"

        function getBlock(arcs::Arcs)::Vi
            block::Vi = Vi()

            curr::Int, next::Int = first(arcs) 

            push!(block, curr)
            deleteat!(arcs, 1)

            while true
                idx = findfirst((i, j)::Pair{Int, Int} -> next == i, arcs)
                idx == nothing && break

                curr, next = arcs[idx]

                push!(block, curr)
                deleteat!(arcs, idx)
            end

            return block
        end

        @debug "Reading blocks profits"
        for (block_id::Int, arcs::Arcs) in block_arcs

            block::Vi = getBlock(arcs)
            push!(data.B, block)

            data.profits[block] = block_profit[block_id]
        end

    end

    # update arcs ids
    data.D.A = collect(keys(data.D.distance))

    # add dummy depot
    @debug "Finding dummy depot"

    depot::Int = data.depot = maximum(collect(keys(data.D.V))) + 1
    data.D.V[depot] = Vertex(depot, -1, -1)

    # dummy weights
    addDummyArcs(data)

    # check feasibility
    @debug "Checking instance connectivity"

    Vb::Si = getBlocksNodes(data)
    distances::ArcCostMap = calculateShortestPaths(data)

    if any((i, j)::Arc -> !in((i, j), keys(distances)), χ(Vb))
        throw(ArgumentError("The SBRP instance it is not connected"))
    end

    # compact in complete graph
    data′::SBRPData, paths′::Dict{Arc, Vi} = createCompleteDigraph(data)

    # check feasibility
    checkArcsFeasibility(data, data′)

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
                                     [(data.depot, i) for i in Vb], # depot arcs
                                     [(i, data.depot) for i in Vb] # depot arcs
                                    )))
    data‴.D.distance = ArcCostMap(map(a::Arc -> a => data.D.distance[a], filter((i, j)::Arc -> i != j, data‴.D.A)))

    # check feasibility
    checkArcsFeasibility(data, data‴)

    # set
    data = data‴

    # compact algorithm
    #  data″, paths″ = create_no_one_degree_paths_digraph(data)

    # check feasibility
    #  check_feasibility(data, data″)

    # return
    #  return data, data′, paths′, data″, paths″
#    return data, data′, paths′
    return data
end

#=
Create a SBRP instance from a Matheus's instance
input:
- app::Dict{String, Any} is the paramenters relation
output:
- data::SBRPData is the SBRP instance
=#
function readSBRPData(app::Dict{String, Any})::SBRPData

    @debug "Reading Matheus's instance"

    depot::Int = 1
    data::SBRPData = SBRPData(
                              InputDigraph(
                                           Dict{Int, Vertex}(1 => Vertex(-1, -1, -1)), 
                                           Arcs(), 
                                           ArcCostMap()
                                          ), 
                              depot,
                              VVi(),
                              parse(Float64, app["vehicle-time-limit"]),
                              ArcCostMap()
                             )
    blocks::Dict{Int, Arcs} = Dict{Int, Arcs}()

    open(app["instance"]) do f::IOStream

        # ignore the first two lines
        readline(f)
        readline(f)

        nNodes::Int = parse(Int, split(readline(f), [' ']; limit=0, keepempty=false)[end])
        nArcs::Int = parse(Int, split(readline(f), [' ']; limit=0, keepempty=false)[end]) + parse(Int, split(readline(f), [' ']; limit=0, keepempty=false)[end])
        nBlocks::Int = parse(Int, split(readline(f), [' ']; limit=0, keepempty=false)[end])

        readline(f)
        readline(f)

        # get arcs
        @debug "Reading arcs"

        for _ in 1:nArcs

            parts::Vector{String} = split(readline(f), ['(', ')', ' ', ',']; limit=0, keepempty=false)

            data.D.distance[Arc(parse(Int, parts[1]), parse(Int, parts[2]))] = floor(Int, parse(Float64, parts[end]))
        end

        # get nodes
        @debug "Reading nodes"

        readline(f)
        for _ in 1:nNodes
            parts::Vector{String} = split(readline(f), [' ']; limit=0, keepempty=false)
            id::Int = parse(Int, parts[1])
            data.D.V[id] = Vertex(id, parse(Float64, parts[2]), parse(Float64, parts[3]))
        end

        # get blocks
        @debug "Reading blocks"

        readline(f)
        for _ in 1:nBlocks
            parts::Vector{String} = split(readline(f), [',', ' ']; limit=0, keepempty=false)
            block::Vi = Vi([parse(Int, part) for part in parts[begin:end - 1]])
            push!(data.B, block)

            # define profit
            data.profits[block] = app["unitary-profits"] ? 1.0 : parse(Float64, parts[end])

        end

    end

    # update arcs ids
    data.D.A = collect(keys(data.D.distance))
#    arcs_set::ArcsSet = ArcsSet(data.D.A)

    # dummy weights
    addDummyArcs(data)

    # check feasibility
    @debug "Checking connectivity"

    Vb::Si = getBlocksNodes(data)

    distances::ArcCostMap = calculateShortestPaths(data)

    if any(a::Arc -> !in(a, keys(distances)), χ(Vb)) 
        throw(InvalidStateException("The SBRP instance it is not connected"))
    end

    # compact in complete graph
    data′::SBRPData, paths′::Dict{Arc, Vi} = createCompleteDigraph(data)

    # check feasibility
    checkArcsFeasibility(data, data′)

#    println(collect(filter(i::Int -> isempty(δ⁻(data′.D.A, i)) || isempty(δ⁺(data′.D.A, i)), collect(keys(data′.D.V)))))
#    println(collect(filter(i::Int -> isempty(δ⁻(data‴.D.A, i)) || isempty(δ⁺(data‴.D.A, i)), collect(keys(data‴.D.V)))))

#    for (a, path) in paths′
#        println("Arc $a")
#        if isempty(path)
#            if !in(a, data.D.A)
#                println("\tArc $a does not exists")
#            end
#        else
#            for arc in zip(path[begin:end - 1], path[begin + 1:end])
#                if !in(arc, data.D.A)
#                    println("\tArc $arc does not exists")
#                end
#            end
#        end
#    end
#    exit(0)

    # consider only shortest paths arcs
    # copy
    data‴ = deepcopy(data)

    # set arcs
    data‴.D.A = collect(ArcsSet(vcat(
                                     [Arc(i, j) for (a::Arc, path::Vi) in paths′ for (i::Int, j::Int) in zip(path[1:end - 1], path[2:end])], # min paths arcs
                                     [Arc(a.first, path[begin]) for (a::Arc, path::Vi) in paths′ if !isempty(path)], # min paths arcs
                                     [Arc(path[end], a.second) for (a::Arc, path::Vi) in paths′ if !isempty(path)], # min paths arcs
                                     [a for (a::Arc, path::Vi) in paths′ if isempty(path)], # min paths arcs (edge case)
                                     [Arc(i, j) for block::Vi in data.B for (i::Int, j::Int) in zip(block[1:end - 1], block[2:end])], # blocks arcs
                                     map(block::Vi -> Arc(block[end], block[begin]), data.B), # blocks arcs
                                     map(i::Int -> Arc(depot, i), collect(Vb)), # depot arcs
                                     map(i::Int -> Arc(i, depot), collect(Vb))  # depot arcs

#                                   [Arc(path[i], path[i + 1]) for (a, path) in paths′ for i in 1:(length(path) - 1)], # min paths arcs
#                                   [Arc(a[1], path[begin]) for (a, path) in paths′ if !∅(path)], # min paths arcs
#                                   [Arc(path[end], a[2]) for (a, path) in paths′ if !∅(path)], # min paths arcs
#                                   [a for (a, path) in paths′ if ∅(path)], # min paths arcs (edge case)
#                                   [Arc(b[i], b[i + 1]) for b in data.B for i in 1:(length(b) - 1)], # blocks arcs
#                                   [Arc(b[end], b[begin]) for b in data.B], # blocks arcs
#                                   [Arc(depot, i) for i in Vb], # depot arcs
#                                   [Arc(i, depot) for i in Vb] # depot arcs

                                    )))
    data‴.D.distance = ArcCostMap(map(a::Arc -> a => data.D.distance[a], data‴.D.A))

    # check feasibility
    checkArcsFeasibility(data, data‴)

    # set
    data = data‴


    # compact algorithm
    #  data″, paths″ = create_no_one_degree_paths_digraph(data)

    # check feasibility
    #  check_feasibility(data, data″)

    # return
    #  return data, data′, paths′, data″, paths″
#    return data, data′, paths′
    return data′
end

#=
Check arcs match of a SBRP instance: This functions checks if the shortest paths of both instances match
input:
- data::SBRPData is the SBRP instance used for basis
- data´::SBRPData is the SBRP that will be verified
output:
- None
exception:
- Case some irregularity is found
=#
function checkArcsFeasibility(data::SBRPData, data′::SBRPData) 

    Vb::Si = getBlocksNodes(data)

    @debug "Calculating shortest paths for the original instance"
    distances::ArcCostMap = calculateShortestPaths(data)

    @debug "Calculating shortest paths for the new instance"
    distances′::ArcCostMap = calculateShortestPaths(data′)

    if any(a::Arc -> !in(a, keys(distances′)), χ(Vb))
        throw(InvalidStateException("The SBRP instance it is not connected", :closed))
    end

    for a::Arc in χ(Vb)
        if abs(distances[a] - distances′[a]) > 1e-3
            throw(InvalidStateException("Original digraph distance of arc $a is $(distances[a]) while the new digraph distance is $(distances′[a])", :closed))
        end
    end
end

#=
Calculate all pair of shortest paths
input:
- data::SBRPData is the SBRP instance
output:
- distances::ArcCostMap is the distances relationship
=#
function calculateShortestPaths(data::SBRPData)::ArcCostMap

    # params
    V::Dict{Int, Vertex} = data.D.V
    A::Arcs = data.D.A
    Vb::Si = getBlocksNodes(data)

    # building adjacency list
    @debug "Building adjacency list"

    adjList::Dict{Int, Vi} = Dict{Int, Vi}(map(i::Int -> i => Vi(), collect(keys(V))))

    for (i::Int, j::Int) in A
        push!(adjList[i], j)
    end

    distances::ArcCostMap = ArcCostMap(map(
                                           a::Arc -> a => typemax(Float64), 
                                           reduce(vcat, map(i::Int -> map(j::Int -> Arc(i, j), collect(keys(V))), collect(Vb)))
                                          ))


    # subrotine to find all shortest paths with root at node i
    function findPaths(i::Int)

        @debug ">> Finding paths from node $i"

        # bfs
        q::Vi = Vi([i])

        distances[Arc(i, i)] = 0.0

        while !isempty(q)

            # curr node
            curr::Int = popfirst!(q)

            # next node
            for next::Int in adjList[curr]

                # cost
                b::Float64 = distances[Arc(i, curr)] + data.D.distance[Arc(curr, next)]

                # update distance and push new nodes
                if !in(next, [i, data.depot]) && distances[Arc(i, next)] > b

                    distances[Arc(i, next)] = b
                    push!(q, next)

                end

            end

        end

    end

    # get paths
    @debug "Obtaining shortest paths in parallel"

    path_tasks::Vector{Task} = Vector{Task}()

    for i::Int in Vb
        t::Task = @task findPaths(i)
        schedule(t)
        push!(path_tasks, t)
    end


    for t::Task in path_tasks
        wait(t)
    end

    return distances
end

#=
Check feasibility of a SBRP complete instance
input:
- data::SBRPData is the SBRP instance
- Vb::Set{Int} is the set of nodes belonging to a block
output:
- None
exception:
- Case some irregularity is found
=#
function checkSBRPCompleteFeasibility(data::SBRPData, Vb::Si)

    @debug "Checking feasibility of the SBRP instance on complete digraph"

    if any(a::Arc -> !in(a, data.D.A), filter((i, j)::Arc -> i != j, χ(Vb)))
        throw(InvalidStateException("The Complete SBRP instance is not feasible"))
    end
end

                                           
#=
Check feasibility of a SBRP instance
input:
- data::SBRPData is the SBRP instance
output:
- None
exception:
- Case some irregulatiry in data´ is found
=#
function checkSBRPfeasibility(data::SBRPData)

    # params
    V′::Si = getBlocksNodes(data)
    A::Arc = data.D.A
    adjList::Dict{Int, Vi} = Dict{Int, Vi}(map(i::Int -> i => Vi(), V′))
    revAdjList::Dict{Int, Vi} = Dict{Int, Vi}(map(i::Int -> i => Vi(), V′))

    # building adjacency list
    @debug "Building adjacency list"

    for (i::Int, j::Int) in A

        push!(adjList[i], j)
        push!(revAdjList[j], i)

    end

    # Kosaraju’s algorithm
    @debug "Kosaraju’s algorithm"

    # outgoing arcs
    @debug "Outgoing arc"

    # bfs
    v::Int = first(V′)
    mandatory_visited_vertices::Si = Si([v])
    visited_vertices::Si = Si([v])
    q::Vi = Vi([v])

    while !isempty(q)

        curr::Int = popfirst!(q) 

        for j::Int in adjList[curr]

            if !in(j, visited_vertices)

                push!(q, j)
                push!(visited_vertices, j)

                if j::Int in V′

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
        for u::Int in V′
            if !in(u, mandatory_visited_vertices)
                throw(InvalidStateException("There is not path between $v and $u"))
            end
        end
    end

    # ingoing arcs
    @debug "Ingoing arc"

    # bfs
    mandatory_visited_vertices = Si([v])
    visited_vertices = Si([v])
    q = [v]
    while !isempty(q)

        curr::Int = popfirst!(q) 

        for i::Int in revAdjList[curr]

            if !in(i, visited_vertices)

                push!(q, i)
                push!(visited_vertices, i)

                if i::Int in V′

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
        for u::Int in V′
            if !in(u, mandatory_visited_vertices)
                throw(InvalidStateException("There is not path between $v and $u"))
            end
        end
    end
end
