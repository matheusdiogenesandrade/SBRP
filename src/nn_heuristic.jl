const getNextBlock(i::Int, data::SBRPData, visited_blocks::VVi)::Tuple{Int, Vi} = argmin(
                                                                                         tuple::Tuple{Float64, Int, Vi} -> first(tuple),
                                                                                         vcat(
                                                                                              #[(time(data, Arc(i, j)), j, block) for block::Vi in data.B for j::Int in block if !in(block, visited_blocks)], 
                                                                                              Vector{Tuple{Float64, Int, Vi}}(reduce(vcat,
                                                                                                                                     map(
                                                                                                                                         block::Vi -> map(j::Int -> Tuple{Float64, Int, Vi}((time(data, Arc(i, j)), j, block)), block), 
                                                                                                                                         setdiff(data.B, visited_blocks)),
                                                                                                                                     init = Vector{Tuple{Float64, Int, Vi}}()
                                                                                                                                    )),
                                                                                              Vector{Tuple{Float64, Int, Vi}}([(typemax(Float64), -1, Vi())]),
                                                                                              Vector{Tuple{Float64, Int, Vi}}([(typemax(Float64), -1, Vi())])
                                                                                             ))[2:end]


function solveNearestNeighborhood(data::SBRPData)::Tuple{Float64, VVi, Vi}
  # setup
  depot::Int = data.depot
  visited_blocks::VVi = VVi()

  # route
  curr::Int = depot
  tour::Vi = Vi([depot])
  incurr_time::Float64 = 0.0
  profit::Float64 = 0.0

  # greedy NN
  while true

    # get next node and block
    next::Int, block::Vi = getNextBlock(curr, data, visited_blocks)                             
    arc::Arc = Arc(curr, next)

    # base case (feasibility check)
    if next == -1 || incurr_time + time(data, arc) + blockTime(data, block) > data.T
        break 
    end

    # save next  
    next != last(tour) && push!(tour, next)

    # mark block as visited
    push!(visited_blocks, block)                                                           

    # update curr and time
    curr = next
    incurr_time += time(data, arc) + blockTime(data, block) 
    profit = profit + data.profits[block] 
  end

  push!(tour, depot)

  return profit, visited_blocks, tour
end
