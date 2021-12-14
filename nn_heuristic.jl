module NearestNeighborhoodHeuristic

using ..Data
using ..Data.SBRP

get_next_block(i::Int64, data::SBRPData, visited_blocks::Set{Array{Int64, 1}}) = min([(SBRP.time(data, (i, j)), j, block) for block in data.B for j in block if !in(block, visited_blocks)]..., (typemax(Float64), -1, []), (typemax(Float64), -1, []))[2:end]


function solve(data::SBRPData)
  # setup
  depot, visited_blocks = data.depot, Set{Array{Int64, 1}}()
  # greedy NN
  curr, tour, time, profit = depot, [depot], 0.0, 0.0
  while true
    (next, block) = get_next_block(curr, data, visited_blocks)                             # get next node and block
    (next == -1 || time + SBRP.time(data, (curr, next)) + time_block(data, block) > data.T) && break # base case (feasibility check)
    push!(tour, next)                                                                      # save next  
    push!(visited_blocks, block)                                                           # mark block as visited
    curr, time, profit = next, time + Data.SBRP.time(data, (curr, next)) + time_block(data, block), profit + data.profits[block] # update curr and time
  end
  return profit, push!(tour, depot)
end

end
