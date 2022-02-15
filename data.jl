module Data

include("symbols.jl")

export Vertex, InpuDigraph, δ⁺, δ⁻, d⁺, d⁻, v⁺, v⁻, time, SBRPData, NORMAL_SPEED

import Unicode

mutable struct Vertex
   id_vertex::Int
   pos_x::Float64
   pos_y::Float64
end

mutable struct InputDigraph # Directed graph
#   V::Vector{Vertex} # set of vertices
   V::Dict{Int, Vertex} # set of vertices
   A::Arcs # set of arcs
   distance::ArcCostMap # time in 40 KM/h
end

δ⁺(A::Arcs, i::Int) = [(j, k) for (j, k) in A if j == i]
δ⁺(A::Arcs, S::Si) = [(i, j) for (i, j) in A if i in S && !in(j, S)]
δ⁺(A::Arcs, S::Vi) = [(i, j) for (i, j) in A if i in S && !in(j, S)]

δ⁺(A::ArcsSet, i::Int) = [(j, k) for (j, k) in A if j == i]
δ⁺(A::ArcsSet, S::Si) = [(i, j) for (i, j) in A if i in S && !in(j, S)]
δ⁺(A::ArcsSet, S::Vi) = [(i, j) for (i, j) in A if i in S && !in(j, S)]

δ⁻(A::Arcs, i::Int) = [(j, k) for (j, k) in A if k == i]
δ⁻(A::Arcs, S::Si) = [(i, j) for (i, j) in A if j in S && !in(i, S)]
δ⁻(A::Arcs, S::Vi) = [(i, j) for (i, j) in A if j in S && !in(i, S)]

δ⁻(A::ArcsSet, i::Int) = [(j, k) for (j, k) in A if k == i]
δ⁻(A::ArcsSet, S::Si) = [(i, j) for (i, j) in A if j in S && !in(i, S)]
δ⁻(A::ArcsSet, S::Vi) = [(i, j) for (i, j) in A if j in S && !in(i, S)]

d⁺(A::ArcsSet, i::Int) = length(δ⁺(A, i))
d⁻(A::ArcsSet, i::Int) = length(δ⁺(A, i))
d⁺(A::Arcs, i::Int) = length(δ⁺(A, i))
d⁻(A::Arcs, i::Int) = length(δ⁺(A, i))

v⁺(A::ArcsSet, i::Int) = [j for (i, j) in δ⁺(A, i)]
v⁻(A::ArcsSet, i::Int) = [j for (j, i) in δ⁻(A, i)]
v⁺(A::Arcs, i::Int) = [j for (i, j) in δ⁺(A, i)]
v⁻(A::Arcs, i::Int) = [j for (j, i) in δ⁻(A, i)]

include("sbrp.jl")

end
