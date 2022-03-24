module Model

include("symbols.jl")

using ..Data
using ..Data.SBRP
using ..NearestNeighborhoodHeuristic
using CPLEX
using JuMP
using BenchmarkTools

export EPS,                   # decimal tolerance rate
       blocks,                # returns all the blocks with at least one node in S
       bfs_sets,              # ...
       add_cuts,              # adds generic cuts (lhs >= rhs) to a given model
       blocks_nodes           # get nodes of a block set

EPS = 1e-4

blocks(B::VVi, S::Si) = VVi([b for b in B if any([i in S for i in b])])
blocks_nodes(B::VVi) = Vi([i for b in B for i in b])

add_cuts(model, cuts) = [@constraint(model, lhs >= rhs) for (lhs, rhs) in cuts]

function bfs(A, i)
  S, q = Set{Int}([i]), [i]
  δ⁺′(j) = [a for a in A if a[1] == j && !in(a[2], S)]
  while !isempty(q)
    curr = popfirst!(q)
    A′ = δ⁺′(curr)
    [(push!(S, next), push!(q, next)) for (curr, next) in A′]
  end
  return S
end

function bfs_sets(A, nodes, source::Int)
  sets = Set{Set{Int}}()
  for i in nodes
    S = bfs(A, i)
    (!in(source, S) && length(S) > 1) && push!(sets, S)
  end
  return sets
end

include("model_sbrp_max.jl")

include("model_sbrp_max_complete.jl")

end
