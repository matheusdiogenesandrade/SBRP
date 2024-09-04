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
