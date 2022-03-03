module Solution

include("symbols.jl")

using ..Data
using ..Data.SBRP
using CPLEX
using JuMP
using DataStructures

export writesol, gettour, check_sbrp_sol, check_atsp_sol, get_info, writeGPX, get_blocks, add_blocks

get_blocks(data::SBRPData, y) = [block for block in data.B if value(y[block]) > 0.5]

get_info(model, data::SBRPData, tour::Vi, B::Vector{Vi}) = Dict{String, String}(
                                       "cost"         => string(objective_value(model)),
                                       "solverTime"   => string(solve_time(model)),
                                       "relativeGAP"  => string(relative_gap(model)),
                                       "nodeCount"    => string(node_count(model)),
                                       "meters"       => string(tour_distance(data, tour)),
                                       "tourMinutes"  => string(tour_time(data, tour, B)),
                                       "blocksMeters" => string(sum(distance_block(data, block) for block in B))
                                      )

function writeGPX(file_path::String, tour::Vector{Vertex})
  open(file_path, "w") do f
    write(f, """<?xml version="1.0" encoding="UTF-8" standalone="no" ?>
    <gpx xmlns="http://www.topografix.com/GPX/1/1" xmlns:gpxx="http://www.garmin.com/xmlschemas/GpxExtensions/v3" xmlns:gpxtpx="http://www.garmin.com/xmlschemas/TrackPointExtension/v1" creator="Oregon 400t" version="1.1" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.topografix.com/GPX/1/1 http://www.topografix.com/GPX/1/1/gpx.xsd http://www.garmin.com/xmlschemas/GpxExtensions/v3 http://www.garmin.com/xmlschemas/GpxExtensionsv3.xsd http://www.garmin.com/xmlschemas/TrackPointExtension/v1 http://www.garmin.com/xmlschemas/TrackPointExtensionv1.xsd">
    <metadata>
    <link href=\"http://www.garmin.com\">
    <text>Garmin International</text>
    </link>
    <time>2009-10-17T22:58:43Z</time>
    </metadata>
    <trk>
    <name>Example GPX Document</name>
    <trkseg>""")
    for v in tour 
      write(f, """<trkpt lat="$(v.pos_y)" lon="$(v.pos_x)">
            </trkpt>\n""")
    end
    write(f, """</trkseg>
    </trk>
    </gpx>""")
  end
end

function writesol(path::String, tour::Vi)
  open(path, "w") do file
    [write(file, "$(v), ") for v in tour]
  end
end

function gettour(data::SBRPData, x, B::Vector{Vi})
  A, depot, tour = [a for a in data.D.A if value(x[a]) > 0.5], data.depot, []
  V = Set{Int}(vcat([i for (i, j) in A], [j for (i, j) in A]))
  #hierholzer's
  adjList, curr_path = Dict{Int, Vi}(i => [j for (i, j) in δ⁺(A, i) for i in 1:Int(floor(value(x[(i, j)]) + 0.5))] for i in V), Stack{Int}()
  curr_path = Stack{Int}()
  push!(curr_path, depot)
  curr_v = first(curr_path)
  while !isempty(curr_path)
    if !isempty(adjList[curr_v])
      push!(curr_path, curr_v)
      curr_v = pop!(adjList[curr_v])
    else
      push!(tour, curr_v)
      curr_v = pop!(curr_path)
    end
  end
  return add_blocks(Vi(reverse(tour)), B)
#  return Vi(reverse(tour))
end

function add_blocks(tour::Vi, B::VVi)
  # add blocks
  tour′, node_blocks = Vi(), Dict{Int, Set{Vi}}(i => Set{Vi}() for b in B for i in b)
  # edge case
  for block in B
    !any([i in tour for i in block]) && error("Street block $block is not being serviced")
  end
  # populate dict
  [push!(node_blocks[j], b) for b in B for j in b]
  for i in tour
    push!(tour′, i)
    (!in(i, keys(node_blocks)) || isempty(node_blocks[i])) && continue # node has no blocks to serve
    for b in node_blocks[i] # for every block
      index = first(findall(x -> x == i, b))
      append!(tour′, [b[j] for j in (index + 1):length(b)]..., [b[j] for j in 1:index]...)
      [delete!(node_blocks[k], b) for k in b] # remove block from remaining nodes
    end
    delete!(node_blocks, i) # remove node from dict
  end
  return tour′
end
#=
function gettour(V::Vi, A::Vector{Tuple{Int, Int}}, depot::Int, x)
#  V, A, depot, tour = V, [a for a in data.D.A if value(x[a]) > 0.5], data.depot, []
  A′, tour = [a for a in A if value(x[a]) > 0.5], Vi()
  #hierholzer's
  adjList, curr_path = Dict{Int, Vi}(i => [j for (i, j) in δ⁺(A′, i) for i in 1:Int(floor(value(x[(i, j)]) + 0.5))] for i in V), Stack{Int}()
  push!(curr_path, depot)
  curr_v = first(curr_path)
  while !isempty(curr_path)
    if !isempty(adjList[curr_v])
      push!(curr_path, curr_v)
      curr_v = pop!(adjList[curr_v])
    else
      push!(tour, curr_v)
      curr_v = pop!(curr_path)
    end
  end
  return Vi(reverse(tour))
end
=#
function check_sbrp_sol(data::SBRPData, tour::Vi, B::VVi)
  V, A = Set{Int}(), Set{Tuple{Int, Int}}(data.D.A)
  # check arcs
  for i in 1:(length(tour) - 1)
    a = (tour[i], tour[i + 1])
    push!(V, a[1])
    push!(V, a[2])
    !in(a, A) && error("Arc $a does not exists")
  end
  # check blocks
  for b in B
    all(!in(i, V) for i in b) && error("Block $b was not served")
  end
end
#=
function check_atsp_sol(tour::Vi, Vb::Dict{Tuple{Int, Vi}, Int}, Vb′::Dict{Tuple{Int, Vi}, Int})
  Vₘ = merge(
    Dict{Int, Vi}(Vb[(i, b)] => b for (i, b) in keys(Vb)),
    Dict{Int, Vi}(Vb′[(i, b)] => b for (i, b) in keys(Vb′))
  )
  i, n = 1, length(tour)
  while i <= n
    node = tour[i]
    if node in keys(Vₘ)
      b, k = Vₘ[node], i + 1
      while k <= n && tour[k] in keys(Vₘ) && Vₘ[tour[k]] == b
        k = k + 1
      end
      k - i + 1 < length(b) * 2 && return false
      i = k
    else
      i = i + 1
    end
  end
  return true
end
=#
end
