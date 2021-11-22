module Solution

using ..Data
using ..Data.SBRP
using CPLEX
using JuMP
using DataStructures

export writesol, gettour, check_sbrp_sol, check_atsp_sol, get_info, writeGPX

get_info(model) = Dict{String, String}(
                                       "cost"         => string(objective_value(model)),
                                       "solverTime"  => string(solve_time(model)),
                                       "relativeGAP" => string(relative_gap(model)),
                                       "nodeCount"   => string(node_count(model))
                                      )

function writeGPX(file_path::String, tour::Array{Vertex, 1})
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

function writesol(path::String, tour::Array{Int64, 1})
  open(path, "w") do file
    [write(file, "$(v), ") for v in tour]
  end
end

function gettour(data::SBRPData, x)
  A, depot, tour = [a for a in data.D.A if value(x[a]) > 0.5], data.depot, []
  V = Set{Int64}(vcat([i for (i, j) in A], [j for (i, j) in A]))
#  println(A)
  #hierholzer's
  adjList, curr_path = Dict{Int64, Vector{Int64}}(i => [j for (i, j) in δ⁺(A, i) for i in 1:Int(floor(value(x[(i, j)]) + 0.5))] for i in V), Stack{Int}()
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
  return Array{Int64}(reverse(tour))
end

#function gettour(data::SBRPData, x)
function gettour(V::Array{Int64}, A::Array{Tuple{Int64, Int64}}, depot::Int64, x)
#  V, A, depot, tour = V, [a for a in data.D.A if value(x[a]) > 0.5], data.depot, []
  A′, tour = [a for a in A if value(x[a]) > 0.5], Array{Int64, 1}()
  #hierholzer's
  adjList, curr_path = Dict{Int64, Vector{Int64}}(i => [j for (i, j) in δ⁺(A′, i) for i in 1:Int(floor(value(x[(i, j)]) + 0.5))] for i in V), Stack{Int}()
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
  return Array{Int64}(reverse(tour))
end

function check_sbrp_sol(data::SBRPData, tour::Array{Int64, 1})
  V′ = Set{Int64}()
  # check arcs
  for i in 1:(length(tour) - 1)
    a = (tour[i], tour[i + 1])
    push!(V′, a[1])
    push!(V′, a[2])
    !in(a, data.D.A) && error("Arc $a does not exists")
  end
  # check blocks
  for b in data.B
    all(!in(i, V′) for i in b) && error("Block $b was not served")
  end
end

function check_atsp_sol(tour::Array{Int64, 1}, Vb::Dict{Tuple{Int64, Array{Int64, 1}}, Int64}, Vb′::Dict{Tuple{Int64, Array{Int64, 1}}, Int64})
  Vₘ = merge(
    Dict{Int64, Array{Int64, 1}}(Vb[(i, b)] => b for (i, b) in keys(Vb)),
    Dict{Int64, Array{Int64, 1}}(Vb′[(i, b)] => b for (i, b) in keys(Vb′))
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

end
