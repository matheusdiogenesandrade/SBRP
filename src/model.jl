using CPLEX
using JuMP
#using BenchmarkTools

const EPS = 1e-4

#=
Get all blocks with some node inside the component ´S´
input: 
- B::Vector{Vector{Int}} is the set of blocks
- S::Set{Int} is the component of nodes
=# 
function blocks(B::VVi, S::Si)::VVi
    return VVi(filter(block::Vi -> !isempty(intersect(S, block)), B))
end

#=
Add cuts to a model
input: 
- model::Model is a Mathematical Programming model
- cuts::Vector{Pair{AffExpr, AffExpr}} is the list of cuts, where each cut has its LHS and RHS
=# 
function addCuts(model::Model, cuts::Vector{Pair{AffExpr, AffExpr}})
    for (lhs::AffExpr, rhs::AffExpr) in cuts
        @constraint(model, lhs >= rhs)
    end
end

#=
Add cuts to a model
input: 
- model::Model is Mathematical Programming model
- cuts::Vector{Pair{AffExpr, AffExpr}} is the list of cuts, where each cut has its LHS and RHS
=# 
#=
function bfs(A::Arcs, i::Int)
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
=#

include("ip_model.jl")
include("cp_model.jl")
