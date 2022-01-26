module NearestNeighborhoodHeuristic

include("symbols.jl")

using ..Data
using ..Data.SBRP

get_next_block(i::Int, data::SBRPData, visited_blocks::Set{Vi}) = min([(SBRP.time(data, (i, j)), j, block) for block in data.B for j in block if !in(block, visited_blocks)]..., (typemax(Float64), -1, []), (typemax(Float64), -1, []))[2:end]


function solve(data::SBRPData)
  # setup
  depot, visited_blocks = data.depot, Set{Vi}()
  curr, tour, time, profit = depot, [depot], 0.0, 0.0

  # greedy NN
  while true
    # get next node and block
    (next, block) = get_next_block(curr, data, visited_blocks)                             

    # base case (feasibility check)
    (next == -1 || time + SBRP.time(data, (curr, next)) + time_block(data, block) > data.T) && break 

    # save next  
    push!(tour, next)                                                                      

    # mark block as visited
    push!(visited_blocks, block)                                                           

    # update curr and time
    curr, time, profit = next, time + Data.SBRP.time(data, (curr, next)) + time_block(data, block), profit + data.profits[block] 
  end

  return profit, visited_blocks, push!(tour, depot)
end

end
