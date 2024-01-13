using Concorde
using Random
using TSPLIB

function solve_tsp(dist_mtx::Matrix{Int})

    n_nodes::Int = size(dist_mtx, 1)

    n_rows::Int = n_nodes
    rows::Vector{String} = Vector{String}()
    for i in 1:n_rows
        push!(rows, join(dist_mtx[i,:], " "))
    end

    name::String = randstring(10)
    filename::String = name * ".tsp"
    open(filename, "w") do io
        write(io, "NAME: $(name)\n")
        write(io, "TYPE: TSP\n")
        write(io, "COMMENT: $(name)\n")
        write(io, "DIMENSION: $(n_nodes)\n")
        write(io, "EDGE_WEIGHT_TYPE: EXPLICIT\n")
        write(io, "EDGE_WEIGHT_FORMAT: FULL_MATRIX \n")
        write(io, "EDGE_WEIGHT_SECTION\n")
        for r::String in rows
            write(io, "$r\n")
        end
        write(io, "EOF\n")
    end

    return __solve_tsp__(filename)

end

function read_solution(filepath)
    sol = readlines(filepath)
    n_nodes = sol[1]
    tour = parse.(Int, split(join(sol[2:end]))) .+ 1
    return tour
end

function tour_length(tour, M)
    n_nodes = length(tour)
    len = 0
    for i in 1:n_nodes
        j = i + 1
        if i == n_nodes
            j = 1
        end

        len += M[tour[i], tour[j]]
    end
    return len
end

function cleanup(name)
    exts = ["mas", "pul", "sav", "sol", "tsp", "res"]
    for ext in exts
        file =  "$(name).$(ext)"
        rm(file, force=true)
        file =  "O$(name).$(ext)"
        rm(file, force=true)
    end
end

function __solve_tsp__(tsp_file::String)

    if !haskey(ENV, "CONCORDE_EXECUTABLE")
        error("CONCORDE_EXECUTABLE variable not set")
    end


    exe = get(ENV, "CONCORDE_EXECUTABLE", nothing)::String
    status = Base.run(`$(exe) $(tsp_file)`, wait=false)

    while !success(status)
        # 
    end    

    name = splitext(basename(tsp_file))[1] :: String
    sol_filepath =  name * ".sol"
    opt_tour = read_solution(sol_filepath) :: Vector{Int}

    tsp = readTSP(tsp_file) :: TSPLIB.TSP
    M = Int.(tsp.weights) :: Matrix{Int}
    opt_len = tour_length(opt_tour, M) :: Int

    cleanup(name)

    return opt_tour, opt_len    
end

function solveTSP(data::SBRPData)::Tuple{SBRPSolution, Dict{String, String}}

    @debug "Solving CSP as TSP"

    # data
    V::Vi = collect(keys(data.D.V))
    B::VVi = data.B

    # convert to TSP
    tsp_distance::Matrix{Int}, 
    V₁_in_relation::Dict{Tuple{Int, Vi}, Int},
    V₁_out_relation::Dict{Tuple{Int, Vi}, Int},
    depot_idx::Int = createTSPInstanceFromCSPInstance(data)

    n_tsp_distance::Int = size(tsp_distance, 1)

    # output
    tour::Vi = Vi()
    len::Int = 0
    info::Dict{String, String} = Dict{String, String}()


    V₁_in_reverse_relation::Dict{Int, Tuple{Int, Vi}} = Dict{Int, Tuple{Int, Vi}}()
    V₁_out_reverse_relation::Dict{Int, Tuple{Int, Vi}} = Dict{Int, Tuple{Int, Vi}}()
    for (i, block) in keys(V₁_in_relation)
        V₁_in_reverse_relation[V₁_in_relation[(i, block)]] = (i, block)
    end
    for (i, block) in keys(V₁_out_relation)
        V₁_out_reverse_relation[V₁_out_relation[(i, block)]] = (i, block)
    end

#    for i::Int in 1:size(tsp_distance, 1)
#        if haskey(V₁_in_reverse_relation, i)
#            println("In: ", V₁_in_reverse_relation[i])
#        elseif haskey(V₁_out_reverse_relation, i)
#            println("Out: ", V₁_out_reverse_relation[i])
#        else
#            println(i, ": error")
#        end
#    end

    # check symmetry
    if tsp_distance != transpose(tsp_distance)

        # data
        T::Int = ceil(data.T)

        #
        sym_tsp_distance::Matrix{Int} = convertATSPToTSP(tsp_distance, T)
        n::Int                        = size(sym_tsp_distance, 1)

        #
        info["solverTime"] = string(@elapsed tour_sym, len = solve_tsp(sym_tsp_distance))
        
#        println(n_tsp_distance)
#        println(tour_sym[iseven.(eachindex(tour_sym))])
#        println(tour_sym[isodd.(eachindex(tour_sym))])
        # check feasibility
        if !all(i::Int -> i > n_tsp_distance, tour_sym[iseven.(eachindex(tour_sym))]) ||
           !all(i::Int -> i <= n_tsp_distance, tour_sym[isodd.(eachindex(tour_sym))])

            error("The tour does not satisfy the symmetric TSP transformation")

        end

        tour, len = tour_sym[isodd.(eachindex(tour_sym))], Int(len + T * n / 2)

        block_idxs = Dict{Vi, Int}(map(tuple -> tuple[2] => tuple[1], enumerate(B)))
        for id::Int in tour
            if id in keys(V₁_in_reverse_relation)
                i, block = V₁_in_reverse_relation[id]
                println(id, ": ", i, " ", block_idxs[block], "(Input)")
            elseif id in keys(V₁_out_reverse_relation)
                i, block = V₁_out_reverse_relation[id]
                println(id, ": ", i, " ", block_idxs[block], "(Output)")
            end
        end

    else

        #
        info["solverTime"] = string(@elapsed tour, len = solve_tsp(tsp_distance))

    end

    # return
    info["cost"] = string(len)

    tour = convertTSPTourToCSPTour(tour, V₁_in_relation, V₁_out_relation, data.depot, depot_idx)
    
    return SBRPSolution(tour, data.B), info

end

function convertTSPTourToCSPTour(
        tsp_tour::Vi, 
        V₁_in_relation::Dict{Tuple{Int, Vi}, Int}, 
        V₁_out_relation::Dict{Tuple{Int, Vi}, Int},
        depot::Int,
        depot_idx::Int
    )::Vi

    #=
    # create reverse relation
    V₁_reverse_relation::Dict{Int, Int} = Dict{Int, Int}()

    for (i::Int, block::Vi) in keys(V₁_in_relation)
        V₁_reverse_relation[V₁_in_relation[(i, block)]] = i
    end

    for (i::Int, block::Vi) in keys(V₁_out_relation)
        V₁_reverse_relation[V₁_out_relation[(i, block)]] = i
    end

    if depot_idx != -1
        V₁_reverse_relation[depot_idx] = depot
    end
    =#

    # output
    tour::Vi = Vi()

    # util
    V₁_in_reverse_relation::Dict{Int, Tuple{Int, Vi}} = Dict{Int, Tuple{Int, Vi}}(map((key, idx)::Pair{Tuple{Int, Vi}, Int} -> idx => key, collect(V₁_in_relation)))
    V₁_out_reverse_relation::Dict{Int, Tuple{Int, Vi}} = Dict{Int, Tuple{Int, Vi}}(map((key, idx)::Pair{Tuple{Int, Vi}, Int} -> idx => key, collect(V₁_out_relation)))

    isInput(idx::Int) = !haskey(V₁_in_reverse_relation, idx) && error("The starting node $idx is not an input node")

    # iterate
    tsp_tour_itr::Iterators.Stateful = Iterators.Stateful(tsp_tour)


    while !isempty(tsp_tour_itr)

        idx::Int = popfirst!(tsp_tour_itr)

        # edge case
        if idx == depot_idx 
            push!(tour, depot)
            continue
        end

        # edge case
        isInput(idx)

        # get node and block
        i::Int, block::Vi = V₁_in_reverse_relation[idx]

        # push
        push!(tour, i)

        # get node block index
        i_block_idx::Int = findfirst(j::Int -> j == i, block)

        # shift block
        shifted_block::Vi = circshift(block, i_block_idx + 1)

        # next
        idx = popfirst!(tsp_tour_itr)

        #
        for j::Int in shifted_block[2:end]
            
            if !haskey(V₁_out_reverse_relation, idx) || V₁_out_reverse_relation[idx] != (j, block) # output case
                error("The node $idx is not the output of node $j of block $shifted_block")
            end

            # next
            idx = popfirst!(tsp_tour_itr)

            if !haskey(V₁_in_reverse_relation, idx) || V₁_in_reverse_relation[idx] != (j, block) # input case
                error("The node $idx is not the input of node $j of block $shifted_block")
            end


            # next
            idx = popfirst!(tsp_tour_itr)
        end

    end

    return unique(tour)
end

function createTSPInstanceFromCSPInstance(data::SBRPData)::Tuple{Matrix{Int}, Dict{Tuple{Int, Vi}, Int}, Dict{Tuple{Int, Vi}, Int}, Int}

    @debug "Creating TSP instance from CSP instance"

    # data
    
    @debug "Get instance data"

    V::Vi                = collect(keys(data.D.V))
    distance::ArcCostMap = data.D.distance
    n::Int               = length(V)
    B::VVi               = data.B
    T::Float64           = data.T

    # new digraph
    
    @debug "Creating TSP instance data"

    V₁_in_relation::Dict{Tuple{Int, Vi}, Int}  = Dict{Tuple{Int, Vi}, Int}() # id => i
    V₁_out_relation::Dict{Tuple{Int, Vi}, Int} = Dict{Tuple{Int, Vi}, Int}() # id => i
    V₁::Vi                                 = Vi() 
    depot_idx::Int                         = -1
    id::Int                                = 1

    # nodes
    
    @debug "Creating TSP nodes"

    for block::Vi in B
        for i::Int in block

            # input
            push!(V₁, id)
            V₁_in_relation[(i, block)]  = id
            id                         += 1

            # output
            push!(V₁, id)
            V₁_out_relation[(i, block)] = id
            id                         += 1

        end
    end

    # depot case
    
    @debug "Depot edge case"

    if data.depot in V
        push!(V₁, id)
        depot_idx = id
    end

    # distances
    
    @debug "Creating TSP distances"

    n₁::Int                   = length(V₁)
    M::Int                    = all(value::Float64 -> value == floor(value), values(distance)) ? 1 : 10000
    all_distances::Int        = M * sum(values(distance)) # an upper bound for the optimal solution suffices
    L::Int                    = all_distances + 1
#    L::Int                    = 0
#    U::Int                    = maximum(map(block::Vi -> length(block), B)) * length(B) * (L + maximum(values(distance))) + 1
    U::Int                    = 1000000
    tsp_distance::Matrix{Int} = fill(U, (n₁, n₁))

    # null cycle case
    
    @debug "Null cycle case"

    for block::Vi in B

        for (i::Int, j::Int) in zip(block[1:end - 1], block[2:end])

            # get nodes
            idx_in_i::Int  = V₁_in_relation[(i, block)] 
            idx_out_j::Int = V₁_out_relation[(j, block)] 

            # set
            tsp_distance[idx_in_i, idx_out_j] = 0

        end

        # closing cycle
        
        for i::Int in block

            # get nodes
            idx_out_i::Int = V₁_out_relation[(i, block)] 
            idx_in_i::Int  = V₁_in_relation[(i, block)] 

            # set
            tsp_distance[idx_out_i, idx_in_i] = 0

        end
       
        idx_in_last::Int   = V₁_in_relation[(last(block), block)] 
        idx_out_first::Int = V₁_out_relation[(first(block), block)] 

        tsp_distance[idx_in_last, idx_out_first] = 0

    end

    #
    
    @debug "Not null distances"

    for block_i::Vi in B
        for block_j::Vi in B

            # edge case
            block_i == block_j && continue

            for i::Int in block_i
                for j::Int in block_j

                    # edge case
#                    i == j && continue

                    # get nodes
                    idx_out_i::Int = V₁_out_relation[(i, block_i)] 
                    idx_in_j::Int  = V₁_in_relation[(j, block_j)] 

                    #
#                    tsp_distance[idx_out_i, idx_in_j] = floor(M * distance[Arc(i, j)]) + U
                    tsp_distance[idx_out_i, idx_in_j] = i == j ? T : T + floor(M * distance[Arc(i, j)])

                end
            end
        end
    end

    # zero diagonal
    for idx::Int in 1:size(tsp_distance, 1)
        tsp_distance[idx, idx] = 0
    end

    # depot case
    
    @debug "Depot edge case"

    if depot_idx != -1

        for block::Vi in B
            for i::Int in block
                # get nodes
                idx_out_i::Int = V₁_out_relation[(i, block)] 
                idx_in_i::Int  = V₁_in_relation[(i, block)] 

                # set
                tsp_distance[depot_idx, idx_in_i]  = floor(M * distance[Arc(data.depot, i)]) + T
                tsp_distance[idx_out_i, depot_idx] = floor(M * distance[Arc(i, data.depot)]) + T
            end
        end

    end

    # return
    return tsp_distance, V₁_in_relation, V₁_out_relation, depot_idx
end

function convertATSPToTSP(tsp_distance::Matrix{Int}, T::Int)::Matrix{Int}

    # data
    n::Int = size(tsp_distance, 1)

    # convert
    sym_distances::Matrix{Int} = fill(1000000, (2 * n, 2 * n))

    for i::Int in 1:n
        for j::Int in 1:n

            value::Int = i == j ? 0 : tsp_distance[i, j]
#            value::Int = i == j ? - M : tsp_distance[i, j]

            sym_distances[j, n + i] = value
            sym_distances[n + i, j] = value

        end
    end

    # return
    return sym_distances

end
