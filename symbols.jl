#=
module Symbols

const ∑ = sum
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
const ≈(a, b) = abs(a - b) <= 1e-3
const ≠(a, b) = a != b 
const ≤(a, b) = a <= b
const ≥(a, b) = a >= b
const ∞ = typemax(Float64)

const flush_println(str) = (println(str); flush(stdout))
const Arc = Tuple{Int, Int}
const Arcs = Vector{Tuple{Int, Int}}
const Vi = Vector{Int}

end
=#

const ∑ = sum
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
const ≈(a, b) = abs(a - b) <= 1e-3
const ≠(a, b) = a != b 
const ≤(a, b) = a <= b
const ≥(a, b) = a >= b
const ∞ = typemax(Float64)

const flush_println(strings...) = (println(strings...); flush(stdout))
const Arc = Tuple{Int, Int}
const Arcs = Vector{Tuple{Int, Int}}
const Vi = Vector{Int}
