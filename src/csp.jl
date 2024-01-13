#=
Create a SBRP instance from an CSP instance
input:
- app::Dict{String, Any} is the paramenters relation
output:
- data::SBRPData is the SBRP instance
=#
function readSBRPDataCSP(app::Dict{String, Any})::SBRPData

    @debug "Reading CSP instance"
    
    depot::Int = -1
    data::SBRPData = SBRPData(
                              InputDigraph(
                                           Dict{Int, Vertex}(), 
                                           Arcs(), 
                                           ArcCostMap()
                                          ), 
                              depot,
                              VVi(),
                              parse(Float64, app["vehicle-time-limit"]),
                              ArcCostMap()
                             )

    open(app["instance"]) do f::IOStream

        # ignore the first line
        readline(f)

        # get number of nodes
        nNodes::Int = parse(Int, split(readline(f), [' ']; limit=0, keepempty=false)[end])

        # store nodes
        for i::Int in 1:nNodes
            data.D.V[i] = Vertex(i, 0.0, 0.0)
        end

        # ignore number of nodes per cluster
        readline(f)

        # ignore distance matrix title
        readline(f)

        # get distances
        @debug "Reading distances"

        for i in 1:nNodes
            distances::Vector{String} = split(readline(f), ['(', ')', ' ', ',']; limit=0, keepempty=false)
            for j in 1:nNodes
                i == j && continue
                data.D.distance[Arc(i, j)] = parse(Int, distances[j])
            end
        end

        # get blocks
        @debug "Reading clusters"

        clusters::VVi = VVi()

        # ignore clusters list title
        readline(f)

        for _ in 1:nNodes

            cluster::Vi = map(k::SubString{String} -> parse(Int, k) + 1, split(readline(f), [',', ' ']; limit=0, keepempty=false))
            push!(clusters, cluster)

        end

        blocks_set::Set{Vi} = Set{Vi}()
        for i::Int in 1:nNodes 

            block::Vi = filter(j::Int -> i in clusters[j], 1:nNodes)

            push!(blocks_set, sort(block))

            # define profit
            data.profits[block] = app["unitary-profits"] ? 1.0 : length(block)

        end

        data.B = collect(blocks_set)

    end

    # update arcs ids
    data.D.A = collect(keys(data.D.distance))

    # check feasibility
    @debug "Checking completeness"

    if any(a::Arc -> !in(a, keys(data.D.distance)), filter((i, j)::Pair{Int, Int} -> i != j, χ(keys(data.D.V))))
        throw("The CSP instance it is not complete ")
    end

    return data
end
