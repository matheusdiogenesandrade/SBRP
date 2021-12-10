module Data

export Vertex, InpuDigraph, δ⁺, δ⁻, time, SBRPData, NORMAL_SPEED
#export Vertex, InpuDigraph, δ⁺, δ⁻, time, time_block, SBRPData, compact, readSBRPData, checkSBRPfeasibility

import Unicode

mutable struct Vertex
   id_vertex::Int64
   pos_x::Float64
   pos_y::Float64
end

mutable struct InputDigraph # Directed graph
   V::Array{Vertex} # set of vertices
   A::Array{Tuple{Int64,Int64}} # set of arcs
   distance::Dict{Tuple{Int64, Int64}, Float64} # time in 40 KM/h
end

δ⁺(A::Array{Tuple{Int64, Int64}}, i::Int64) = [(j, k) for (j, k) in A if j == i]
δ⁺(A::Array{Tuple{Int64, Int64}}, S::Set{Int64}) = [(i, j) for (i, j) in A if i in S && !in(j, S)]
δ⁺(A::Array{Tuple{Int64, Int64}}, S::Array{Int64}) = [(i, j) for (i, j) in A if i in S && !in(j, S)]
δ⁻(A::Array{Tuple{Int64, Int64}}, i::Int64) = [(j, k) for (j, k) in A if k == i]
δ⁻(A::Array{Tuple{Int64, Int64}}, S::Set{Int64}) = [(j, k) for (j, k) in A if !in(k, S)]
δ⁻(A::Array{Tuple{Int64, Int64}}, S::Array{Int64}) = [(j, k) for (j, k) in A if !in(k, S)]
NORMAL_SPEED = (40.0 * 10^3)/60.0
time(data, a) = data.D.distance[a] / NORMAL_SPEED

include("sbrp.jl")

end
