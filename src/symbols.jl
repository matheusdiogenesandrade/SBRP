const ∧(x...) = all(x)
const ∨(x...) = any(x)
const →(a, b) = !a ∨ b
const ←(a, b) = !b ∨ a
const ∩(x...) = intersect(x...)
const ∪(x...) = union(x...)
const ⊆(a, b) = all(x -> x ∈ b, a)
const ⊂(a, b) = a ⊆ b ∧ !all(x -> x ∈ a, b)
const ⊇(a, b) = b ⊆ a
const ⊃(a, b) = b ⊂ a
const ≈(a::Real, b::Real) = abs(a - b) <= 1e-3
const ≠(a, b) = a != b 
const ≤(a::Real, b::Real) = a <= b
const ≥(a::Real, b::Real) = a >= b
const ∞ = typemax(Float64)
const χ(x) = [Pair(i, j) for i in x for j in x]
const χ(x, y) = [Pair(i, j) for i in x for j in y]

const ∅(x) = isempty(x)

const Arc = Pair{Int, Int}
const Arcs = Vector{Arc}
const Arcs(arcs::Vector{Tuple{Int, Int}}) = Arcs(map((i, j)::Tuple{Int, Int} -> (i => j), arcs)) 
const ArcsSet = Set{Arc}
const Vi = Vector{Int}
const VVi = Vector{Vi}
const Si = Set{Int}
const ArcCostMap = Dict{Arc, Float64}
