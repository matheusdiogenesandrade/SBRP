#=
getBlocks(data::SBRPData, y) = [block for block in data.B if value(y[block]) > 0.5]

getInfo(model, data::SBRPData, tour::Vi, B::Vector{Vi}) = Dict{String, String}(
"cost"         => string(has_values(model) ? objective_value(model) : sum(block -> data.profits[block], B)),
"solverTime"   => string(has_values(model) ? solve_time(model) : 3600),
"relativeGAP"  => string(has_values(model) ? relative_gap(model) : "-"),
"nodeCount"    => string(has_values(model) ? node_count(model) : 0),
"meters"       => string(tour_distance(data, tour)),
"tourMinutes"  => string(tour_time(data, tour, B)),
"blocksMeters" => string(sum(distance_block(data, block) for block in B))
)
=#


#=
Write GPX file
input:
- file_path::String is the file path the GPX solution file will be saved
- tour::Vector{Vertex} is the route to be written in a file
output:
- None
=#
#=
function writeGPX(file_path::String, tour::Vector{Vertex})

    @debug "Writing GPX route file"

    open(file_path, "w") do f::IOStream

        write(f, """<?xml version="1.0" encoding="UTF-8" standalone="no" ?>
              <gpx xmlns="http://www.topografix.com/GPX/1/1" xmlns:gpxx="http://www.garmin.com/xmlschemas/GpxExtensions/v3" xmlns:gpxtpx="http://www.garmin.com/xmlschemas/TrackPointExtension/v1" creator="Oregon 400t" version="1.1" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.topografix.com/GPX/1/1 http://www.topografix.com/GPX/1/1/gpx.xsd http://www.garmin.com/xmlschemas/GpxExtensions/v3 http://www.garmin.com/xmlschemas/GpxExtensionsv3.xsd http://www.garmin.com/xmlschemas/TrackPointExtension/v1 http://www.garmin.com/xmlschemas/TrackPointExtensionv1.xsd">
              <metadata>
              <link href=\"http://www.garmin.com\">
              <text>Garmin International</text>
              </link>
              <time>2009-10-17T22:58:43Z</time>
              </metadata>
              <trk>
              <name>Example GPX Document</name>
              <trkseg>""")

        for v::Vertex in tour 

            write(f, "<trkpt lat="$(v.pos_y)" lon="$(v.pos_x)"></trkpt>\n")

        end

        write(f, "</trkseg> </trk> </gpx>")

    end
end
=#

#=
Write solution in a file
input:
- file_path::String is the file path the solution file will be saved
- tour::Vector{Vertex} is the route to be written in a file
- data::SBRPData is the SBRP instance in which the tour was built on
- B::Vector{Vector{Int}} is the set of blocks serviced by the route ´route´
output:
- None
=#
function writeSolution(file_path::String, data::SBRPData, solution::SBRPSolution)

    @debug "Writing solution in a file"

    # get solution params
    tour::Vi = solution.tour

    serviced_blocks::VVi = solution.B

    # write route
    tour_no_depot::Vi = filter(i::Int -> i != data.depot, tour)

    sol_string::String = join(tour_no_depot, ", ") * "\n"

    for block::Vi in serviced_blocks
        sol_string *= join(block, ", ") * "\n"
    end

    write(open(file_path * ".sol", "w"), sol_string)

    # write GPS
#    writeGPX(file_path * ".gpx", map(i::Int -> data.D.V[i], tour_no_depot)) 
end

#=
Check SBRP solution feasibility
input:
- data::Vector{Vertex} is the route to be written in a file
output:
- None
=#
function checkSBRPSolution(data::SBRPData, solution::SBRPSolution)

    @debug "Checking the feasibility of a SBRP solution"

    # instance parameters
    A::ArcsSet = ArcsSet(data.D.A)

    # solution parameters
    tour::Vi = solution.tour
    visited_nodes::Si = Set{Int}(tour)
    serviced_blocks::VVi = solution.B

    # check arcs
    for (i::Int, j::Int) in zip(tour[begin:end - 1], tour[begin + 1:end])

        if !in(Arc(i, j), A) 

            throw(ErrorException("Arc $(Arc(i, j)) does not exists"))
        end

    end

    # check blocks
    for block::Vi in serviced_blocks

        if all(i::Int -> !in(i, visited_nodes), block)
            throw(ErrorException("Block $block was not served"))
        end

    end

end

#=
Reads a SBRP solution from a file
input:
- file_path::String is the file directory the solution file is stored at
- data::SBRPData is the SBRP instance
output:
- (route, B)::Pair{Vector{Int}, Vector{Vector{Int}}} is the pair containing the route and serviced blocks
=#
function readSolution(file_path::String, data::SBRPData)::SBRPSolution
    # get lines
    lines::Vector{String} = readlines(file_path)

    # get route
    tour::Vi = map(node -> parse(Int, node), split(first(lines), ", ", keepempty=false))
    pushfirst!(tour, data.depot)
    push!(tour, data.depot)

    # get blocks
    serviced_blocks::VVi = map(block_line::String -> map(node::Int -> parse(Int, node), split(block_line, ", ", keepempty=false)), lines[2:end])

    # return
    return  SBRPSolution(tour, serviced_blocks)
end


#=
function add_blocks(tour::Vi, B::VVi)
    # add blocks
    tour′, node_blocks = Vi(), Dict{Int, Set{Vi}}(i => Set{Vi}() for b in B for i in b)
    # edge case
    for block in B
        !any([i in tour for i in block]) && error("Street block $block is not being serviced")
    end
    # populate dict
    [push!(node_blocks[j], b) for b in B for j in b]
    for i in tour
        push!(tour′, i)
        (!in(i, keys(node_blocks)) || isempty(node_blocks[i])) && continue # node has no blocks to serve
        for b in node_blocks[i] # for every block
            index = first(findall(x -> x == i, b))
            append!(tour′, [b[j] for j in (index + 1):length(b)]..., [b[j] for j in 1:index]...)
            [delete!(node_blocks[k], b) for k in b] # remove block from remaining nodes
        end
        delete!(node_blocks, i) # remove node from dict
    end
    return tour′
end

function gettour(data::SBRPData, x, B::Vector{Vi})
    A, depot, tour = [a for a in data.D.A if value(x[a]) > 0.5], data.depot, []
    V = Set{Int}(vcat([i for (i, j) in A], [j for (i, j) in A]))
    #hierholzer's
    adjList, curr_path = Dict{Int, Vi}(i => [j for (i, j) in δ⁺(A, i) for i in 1:Int(floor(value(x[(i, j)]) + 0.5))] for i in V), Stack{Int}()
    curr_path = Stack{Int}()
    push!(curr_path, depot)
    curr_v = first(curr_path)
    while !isempty(curr_path)
        if !isempty(adjList[curr_v])
            push!(curr_path, curr_v)
            curr_v = pop!(adjList[curr_v])
        else
            push!(tour, curr_v)
            curr_v = pop!(curr_path)
        end
    end
    return add_blocks(Vi(reverse(tour)), B)
    #  return Vi(reverse(tour))
end
=#
#=
function gettour(V::Vi, A::Vector{Tuple{Int, Int}}, depot::Int, x)
#  V, A, depot, tour = V, [a for a in data.D.A if value(x[a]) > 0.5], data.depot, []
A′, tour = [a for a in A if value(x[a]) > 0.5], Vi()
#hierholzer's
adjList, curr_path = Dict{Int, Vi}(i => [j for (i, j) in δ⁺(A′, i) for i in 1:Int(floor(value(x[(i, j)]) + 0.5))] for i in V), Stack{Int}()
push!(curr_path, depot)
curr_v = first(curr_path)
while !isempty(curr_path)
if !isempty(adjList[curr_v])
push!(curr_path, curr_v)
curr_v = pop!(adjList[curr_v])
else
push!(tour, curr_v)
curr_v = pop!(curr_path)
end
end
return Vi(reverse(tour))
end
=#
#=
function check_atsp_sol(tour::Vi, Vb::Dict{Tuple{Int, Vi}, Int}, Vb′::Dict{Tuple{Int, Vi}, Int})
Vₘ = merge(
Dict{Int, Vi}(Vb[(i, b)] => b for (i, b) in keys(Vb)),
Dict{Int, Vi}(Vb′[(i, b)] => b for (i, b) in keys(Vb′))
)
i, n = 1, length(tour)
while i <= n
node = tour[i]
if node in keys(Vₘ)
b, k = Vₘ[node], i + 1
while k <= n && tour[k] in keys(Vₘ) && Vₘ[tour[k]] == b
k = k + 1
end
k - i + 1 < length(b) * 2 && return false
i = k
else
i = i + 1
end
end
return true
end
=#

#=
Read solution file
input:
- file_path::String is the file path where the solution is located at
- data::SBRPData is the SBRP instance
output:
- (route, B)::Pair{Vi, VVi} is a pair, where ´route´ is the tour, and ´B´ is the list of blocks
=#
