module Data

include("symbols.jl")

export Vertex, InpuDigraph, δ⁺, δ⁻, time, SBRPData, NORMAL_SPEED

import Unicode

mutable struct Vertex
   id_vertex::Int
   pos_x::Float64
   pos_y::Float64
end

mutable struct InputDigraph # Directed graph
   V::Vector{Vertex} # set of vertices
   A::Vector{Arc} # set of arcs
   distance::Dict{Arc, Float64} # time in 40 KM/h
end

δ⁺(A::Arcs, i::Int) = [(j, k) for (j, k) in A if j == i]
δ⁺(A::Arcs, S::Set{Int}) = [(i, j) for (i, j) in A if i in S && !in(j, S)]
δ⁺(A::Arcs, S::Vi) = [(i, j) for (i, j) in A if i in S && !in(j, S)]
δ⁻(A::Arcs, i::Int) = [(j, k) for (j, k) in A if k == i]
δ⁻(A::Arcs, S::Set{Int}) = [(j, k) for (j, k) in A if !in(k, S)]
δ⁻(A::Arcs, S::Vi) = [(j, k) for (j, k) in A if !in(k, S)]

include("sbrp.jl")

end
