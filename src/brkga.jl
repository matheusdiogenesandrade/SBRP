using BrkgaMpIpr
using ConfParser
using Dates
using Printf

"""
    @enum StopRule

Controls stop criteria. Stops either when:
- a given number of `GENERATIONS` is given;
- or a `TARGET` value is found;
- or no `IMPROVEMENT` is found in a given number of iterations.
"""
@enum StopRule begin
    GENERATIONS = 0
    TARGET = 1
    IMPROVEMENT = 2
end

"""
    function parse(::Type{StopRule}, value::String)::StopRule

Parse `value` into a `StopRule`.
"""
function parseRule(::Type{StopRule}, value::String)::StopRule

    @debug "Parsing GA rules"

    local_value::Char = uppercase(strip(value)[1])

    dict::Dict{Char, StopRule} = Dict{Char, StopRule}('G' => GENERATIONS, 'T' => TARGET, 'I' => IMPROVEMENT)

    local_value in keys(dict) && return dict[local_value]

    throw(ArgumentError("cannot parse $value as StopRule"))
end

function dijkstra(data::SBRPData, idxs_blocks::Vi)::Tuple{Float64, Dict{Tuple{Int, Int}, Float64}}

    @debug "Dijkstra algorithm"


    # attrs
    
    @debug "Getting params"

    m::Int = length(idxs_blocks)
    B::VVi         = data.B
    A::Arcs        = data.D.A
    T::Float64     = data.T
    first_idx::Int = first(idxs_blocks) 

    # dijkstra distance <<idx_block, node>, min path time> stores the min path time from the depot to the node `node' in the block at `idx_block'
    consumed_times::Dict{Tuple{Int, Int}, Float64} = Dict{Tuple{Int, Int}, Float64}(reduce(vcat,
                                                                                           map(idx_block::Int -> 
                                                                                               map(i::Int -> (idx_block, i)::Tuple{Int, Int} => T + 1, B[idx_block]),
                                                                                               idxs_blocks
                                                                                              )
                                                                                          ))

    # initial times

    @debug "Initial times"

    first_block_time::Float64 = blockTime(data, B[first_idx])
    for i::Int in B[first_idx]
        consumed_times[(first_idx, i)] = first_block_time
    end

    # BFS

    @debug "BFS"

    for i::Int in 1:m - 1

        # get indexes

        @debug "Get indexes"

        idx_curr_block::Int = idxs_blocks[i]
        idx_next_block::Int = idxs_blocks[i + 1] 

        # get curr and next blocks

        @debug "Get curr and next block"

        curr_block::Vi = B[idx_curr_block]
        next_block::Vi = B[idx_next_block]

        # calculate next block service time
        
        @debug "Get next block service time"

        next_block_time::Float64 = blockTime(data, next_block)

        # intersecting nodes
       
        @debug "Get intersecting nodes"

        intersection::Si = Si(intersect(curr_block, next_block))

        for j::Int in intersection
            consumed_times[(idx_next_block, j)] = next_block_time + consumed_times[(idx_curr_block, j)]
        end

        # non intersecting nodes
        
        @debug "Get non-intersecting nodes"

        for i::Int in curr_block
            for j::Int in next_block
                if Arc(i, j) in keys(data.D.distance)

                    consumed_times[(idx_next_block, j)] = min(consumed_times[(idx_next_block, j)], consumed_times[(idx_curr_block, i)] + time(data, Arc(i, j)) + next_block_time)

                end
            end
        end
    end

    @debug "Get last block index"

    last_block_idx::Int = idxs_blocks[findlast(idx::Int -> any(i::Int -> consumed_times[(idx, i)] <= T, B[idx]), idxs_blocks)]

    @debug "Return"

    candidates::Vi = filter(i::Int -> consumed_times[(last_block_idx, i)] <= T, B[last_block_idx])

    return minimum(map(i::Int -> consumed_times[(last_block_idx, i)], candidates)), consumed_times
end

function getLongestProfitRoute(data::SBRPData, consumed_times::Dict{Tuple{Int, Int}, Real}, idxs_blocks::Vi)

    #=
    # params
    V = data.D.V
    depot = data.depot
    B = data.B

    # blocks order
    blocks_order::VVi = map(idx::Int -> B[idx], idxs_blocks)
    pushfirst!(blocks_order, [depot])

    # build DAG
    # build nodes
    dag_V::Dict{Int, Vertex} = Dict{Int, Vertex}()

    block_node_to_id::Dict{Tuple{Vi, Int}, Int} = Dict{Tuple{Vi, Int}, Int}()
    id_to_block_node::Dict{Int, Tuple{Vi, Int}} = Dict{Int, Tuple{Vi, Int}}()

    id::Int = 0
    for block in blocks_order
    for i in block
    dat_V[id] = V[i]
    block_node_to_id[(block, i)] = id
    id_to_block_node[id] = (block, i)

    id += 1
    end
    end

    # build arcs
    dag_A::Arcs = Arcs()
    enumeration = enumerate(blocks_order)
    for idx_prev, block_prev in enumeration[1:end - 1]
    for _, block_next in enumeration[idx_prev + 1:end]
    for i in block_prev
    for j in block_next
    node_prev = block_node_to_id[(block_prev, i)]
    node_next = block_node_to_id[(block_next, j)]

    push!(dag_arcs, Arc(node_prev, node_next))
    end
    end
    end
    end

    # build revAdjList
    revAdjList::Dict{Int, Vi} = Dict{Int, Vi}(id => Vi())
    for (id_i, id_j) in dag_A 
    push!(revAdjList[id_j], id_i)
    end

    # times and profits
    times::ArcCostMap = ArcCostMap()
    profits::ArcCostMap = ArcCostMap()
    for (id_i, id_j) in dag_A
    i, block_i = id_to_block_node[id_i]
    j, block_j = id_to_block_node[id_j]

    times[(id_i, id_j)] = time(data, (i, j)) + time_block(data, block_j)
    profits[(id_i, id_j)] = data.profits[block_j]
    end

    # get topological sorting
    topological_order::Vi = Vi()
    for block in blocks_order
    for i in block
    id_i = block_node_to_id[(block, i)]

    push!(topological_order, id_i)
    end
    end

    #
    path::Vi = Vi()
    pred::Dict{Int, Int} = Dict{Int, Int}(i => i for i in keys(dag_V))
    incurr_time::Vector{Float64} = Vector{Float64}(zeros(length(pred)))
    incurr_profit::Vector{Float64} = Vector{Float64}(zeros(length(pred)))

    # get distances
    for v in topological_order
    for u in revAdjList[v]

    newdist = dist[u] + times[(u, v)]
    newprofit = 

    if newdist > dist[v]
    dist[v] = newdist
    pred[v] = u
    end

    end
    end

    # retrieve path 
    v = argmax(x -> dist[x], vertices(g))
    push!(path, v)

    while pred[v] != v
    v = pred[v]
    push!(path, v)
    end

    # return 
    return reverse(path)
    =#

end

function getDijkstraRoute(data::SBRPData, consumed_times::Dict{Tuple{Int, Int}, Float64}, idxs_blocks::Vi)::Vi

    @debug "Get dijkstra route"

    # attrs

    @debug "Get params"

    m::Int = findlast(idx::Int -> any(i::Int -> consumed_times[(idx, i)] <= data.T, data.B[idx]), idxs_blocks)

    # edge case
    m == 0 && return tour

    tour::Vi       = Vi()
    B::VVi         = data.B
    A::Arcs        = data.D.A
    T::Float64     = data.T
    first_idx::Int = first(idxs_blocks)
    last_idx::Int  = idxs_blocks[m]

    # initialize tour
    
    @debug "Initialize tour"

    push!(tour, B[last_idx][findmin(i -> consumed_times[(last_idx, i)], B[last_idx])[2]])

    # BFS
    
    @debug "BFS"

    for i::Int in reverse(2:m)

        # get indexes

        @debug "Get curr and next block indexes"

        idx_curr_block::Int = idxs_blocks[i]
        idx_prev_block::Int = idxs_blocks[i - 1]

        # get curr and prev blocks

        @debug "Get curr and next block"

        curr_block::Vi = B[idx_curr_block]
        prev_block::Vi = B[idx_prev_block]

        # get intersection

        @debug "Get intersection"

        intersection::Vi = ∩(curr_block, prev_block)

        # get block time

        @debug "Get curr block time"

        curr_block_time::Float64 = blockTime(data, B[idx_curr_block])

        # if no intersecting nodes push backwards candidate

        @debug "No intersecting nodes edge-case"

        j::Int = last(tour)
        if !in(j, intersection)
            candidates::Vi = filter(
                                    i::Int -> consumed_times[(idx_prev_block, i)] + time(data, Arc(i, j)) + curr_block_time <= consumed_times[(idx_curr_block, j)], 
                                    prev_block
                                   )
            push!(tour, first(candidates))
        end

    end

    return reverse(tour)

end

"""
Countings
"""
N_FEASIBLE::Int        = 0
N_INFEASIBLE::Int      = 0
DECODING_TIME::Float64 = 0

function getIdxsBlocks(chromosome::Array{Float64}, data::SBRPData)::Vi
    # attrs
    m::Int = length(data.B)

    # get all the alleles with weight > 0.5 and sort them in reversed order 
    permutation::Array{Tuple{Float64, Int64}} = Array{Tuple{Float64, Int64}}(undef, m)

    for (idx::Int, key::Float64) in enumerate(chromosome)
        permutation[idx] = Tuple{Float64, Int}((key, idx))
    end

    sort!(permutation, rev = true)

    # consider only blocks with allele > 0.5
    #  filter!(s -> s[1] > 0.5, permutation)

    # return indexes of the selected blocks
    return map((key, idx)::Tuple{Float64, Int} -> idx, permutation)
end

function decode!(chromosome::Array{Float64}, data::SBRPData, rewrite::Bool)::Float64
    # import blogal variables
    global N_FEASIBLE
    global N_INFEASIBLE 
    global DECODING_TIME 

    DECODING_TIME += curr_decode_time::Float64 = @elapsed begin

        # attrs
        B::VVi = data.B

        # get min tour time
        idxs_blocks::Vi = getIdxsBlocks(chromosome, data)
        tour_time::Float64, consumed_times::Dict{Tuple{Int, Int}, Float64} = dijkstra(data, idxs_blocks)

        idxs_blocks = idxs_blocks[begin:findlast(idx::Int -> any(i::Int -> consumed_times[(idx, i)] <= data.T, data.B[idx]), idxs_blocks)]
    end

    if tour_time > data.T # check feasibility

        throw(InvalidStateException("Not here"))
        #      println("Infeasile with time ", tour_time)
        N_INFEASIBLE += 1
        return -∞

    else # return profit
        #      println("Feasile with time ", tour_time)
        N_FEASIBLE += 1
        return sum(idx::Int -> data.profits[B[idx]], idxs_blocks)

    end
end
