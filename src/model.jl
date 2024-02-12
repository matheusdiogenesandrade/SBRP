using CPLEX
using JuMP
using DataStructures
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

#=
Create dominance digraph
input: 
- data::SBRPData is the set of blocks
output:
- dominance_digraph::InputDigraph is the dominance digraph
=# 
function getDominanceDigraph(data::SBRPData)::InputDigraph

    # data
    V::Dict{Int, Vertex}             = data.D.V 
    B::VVi                           = data.B
    depot::Int                       = data.depot
    node_blocks::Dict{Int, VVi}      = Dict{Int, VVi}(map(i::Int -> i => filter(block::Vi -> i in block, B), collect(keys(V))))
    node_blocks_union::Dict{Int, Vi} = Dict{Int, Vi}(map(i::Int -> i => reduce(union, node_blocks[i]; init=[]), collect(keys(V))))

    # digraphs
    dominance_digraph::InputDigraph = InputDigraph(
                                                   V, 
                                                   Arcs(), 
                                                   ArcCostMap()
                                                  )

    # get dominance relationship (pair-per-pair)
    for i::Int in keys(V)

        # i data
        node_i_blocks::VVi      = node_blocks[i] 
        node_i_blocks_union::Vi = node_blocks_union[i] 

        for j::Int in keys(V)

            # edge case
            (i == j || depot in [i, j]) && continue

            # j data
            node_j_blocks::VVi      = node_blocks[j]
            node_j_blocks_union::Vi = node_blocks_union[j] 

            # does not dominate
            !issubset(node_j_blocks, node_i_blocks) && continue

            #
            push!(dominance_digraph.A, Arc(i, j))
        end
    end

    #
    return dominance_digraph
end

#=
Create a compacted dominance digraph, which is a DAG
input: 
- dominance_digraph::InputDigraph is the dominance digraph
output:
- compact_dominance_digraph::InputDigraph is the DAG dominance digraph
- root::Int is the root node of the resulting DAG
- nodes_sccs::Dict{Int, Vi} is the relation mapping each node of compact_dominance_digraph to a complete subcomponent of dominance_digraph
=# 
function getCompactDominanceDigraph(dominance_digraph::InputDigraph)::Tuple{InputDigraph, Int, Dict{Int, Vi}}

    # data
    V::Dict{Int, Vertex} = dominance_digraph.V
    A::Arcs              = dominance_digraph.A
    arcs_set::ArcsSet    = ArcsSet(A)
    adjList::Dict{Int, Vi} = Dict{Int, Vi}(map(i::Int -> i => Vi(), collect(keys(V))))
    for (i::Int, j::Int) in A
        push!(adjList[i], j)
    end
    SCCs::Set{Si} = Set{Si}()

    # digraph
    compact_dominance_digraph::InputDigraph = InputDigraph(
                                                           Dict{Int, Vertex}(), 
                                                           Arcs(), 
                                                           ArcCostMap()
                                                          )
    nodes_sccs::Dict{Int, Vi} = Dict{Int, Vi}()
    
    # Tarjan's SCC
    visited_nodes::Si      = Si()
    processed::Si          = Si()
    num::Dict{Int, Int}    = Dict{Int, Int}()
    lowest::Dict{Int, Int} = Dict{Int, Int}()
    counter::Int           = 1
    stack::Stack{Int}      = Stack{Int}() 

    function DFS(i::Int) 

        # pre steps
        num[i] = counter
        lowest[i] = num[i]
        counter += 1
        push!(visited_nodes, i)
        push!(stack, i)

        # explore
        for j::Int in adjList[i]

            if !in(j, visited_nodes) # if not visited
               
                DFS(j)
                lowest[i] = min(lowest[i], lowest[j])

            elseif !in(j, processed) # if not processed
                
                lowest[i] = min(lowest[i], num[j])

            end

        end

        push!(processed, i)

        if lowest[i] == num[i]

            scc::Vi = Vi()
            scc_node::Int = pop!(stack)

            while scc_node != i
                push!(scc, scc_node)
                scc_node = pop!(stack)
            end

            push!(scc, scc_node)

            # store
                push!(SCCs, Si(scc))

        end
    end

    #
    for i::Int in keys(V)
        DFS(i)
    end

    # preprocess SCCs
    filter!(scc::Si -> all(scc1::Si -> scc1 == scc || !issubset(scc, scc1), SCCs), SCCs)

    # prepare output
    # nodes
    
    # root
    root::Int = 0
    nodes_sccs[root] = []
    compact_dominance_digraph.V[root] = Vertex(root, 0.0, 0.0)

    # SCCs
    for (idx::Int, scc::Si) in enumerate(SCCs)
        nodes_sccs[idx] = collect(scc)
        compact_dominance_digraph.V[idx] = Vertex(idx, 0.0, 0.0)
    end

    # arcs
    compact_dominance_digraph_nodes::Vi = collect(keys(compact_dominance_digraph.V))
    for (i::Int, j::Int) in χ(compact_dominance_digraph_nodes, compact_dominance_digraph_nodes)
        # edge case
        (i == root || i == j) && continue

        #
        a::Arc = Arc(i, j)
        
        # get SCCs
        scc_i::Vi = nodes_sccs[i]
        scc_j::Vi = nodes_sccs[j]
        
        # cases
        if any((k, l)::Arc -> Arc(k, l) in arcs_set, χ(scc_i, scc_j)) # case arc exists

            push!(compact_dominance_digraph.A, a)
            compact_dominance_digraph.distance[a] = length(setdiff(scc_j, scc_i))

        end
    end

    for a in compact_dominance_digraph.A
        println(a)
    end

    # from root arcs
    # revAdjList
    compact_digraph_revAdjList::Dict{Int, Vi} = Dict{Int, Vi}(map(j::Int -> j => Vi(), compact_dominance_digraph_nodes))
    for (i::Int, j::Int) in compact_dominance_digraph.A
        push!(compact_digraph_revAdjList[j], i)
    end

    # root nodes
    for j::Int in filter(j::Int -> j != root && isempty(compact_digraph_revAdjList[j]), compact_dominance_digraph_nodes)
        a::Arc = Arc(root, j)
        push!(compact_dominance_digraph.A, a)
        compact_dominance_digraph.distance[a] = length(nodes_sccs[j]) 
    end

    return compact_dominance_digraph, root, nodes_sccs

end

function getDAGTransitiveReduction(compact_dominance_digraph::InputDigraph, root::Int)::InputDigraph

    # data
    V::Dict{Int, Vertex} = compact_dominance_digraph.V
    A::Arcs = compact_dominance_digraph.A

    # util
    adjList::Dict{Int, Vi} = Dict{Int, Vi}(map(i::Int -> i => Vi(), collect(keys(V))))

    for (i::Int, j::Int) in A
        push!(adjList[i], j)
    end

    # output
    intransitive_dag::InputDigraph = InputDigraph(V,
                                                  Arcs(),
                                                  ArcCostMap())

    # topological sorting
    topological_sorting::Vi = Vi()
    permanent_mark::Si = Si()
    buffer::Si = Si()

    function DFSForTopologicalSort(i::Int)

        # edge case
        i in permanent_mark && return

        # store in buffer
        push!(buffer, i)

        for j::Int in adjList[i]
            DFSForTopologicalSort(j)
        end

        # remove from buffer
        delete!(buffer, i)

        # store as permanent
        push!(permanent_mark, i)

        # store in list
        push!(topological_sorting, i)

    end

    DFSForTopologicalSort(root)

    reverse!(topological_sorting)

    # transitive reduction
    topological_sorting_idx::Dict{Int, Int} = Dict{Int, Int}(map((idx, i)::Tuple{Int, Int} -> i => idx, enumerate(topological_sorting)))
    intransitive_arcs::Arcs = Arcs()
    reachable_nodes::Dict{Int, Si} = Dict{Int, Si}()

    for i::Int in reverse(topological_sorting)

        # default
        reachable_nodes[i] = Si([i])

        for j::Int in sort(adjList[i]; by = j::Int -> topological_sorting_idx[j])

            # edge case
            in(j, reachable_nodes[i]) && continue

            #
            push!(intransitive_arcs, Arc(i, j))
            union!(reachable_nodes[i], reachable_nodes[j])
        end
    end

    # set
    intransitive_dag.A = intransitive_arcs

    # filter distances
    intransitive_dag.distance = ArcCostMap(map(a::Arc -> a => compact_dominance_digraph.distance[a], intransitive_dag.A))

    # return
    return intransitive_dag
end

function getMaximalPathsDAG(intransitive_dag::InputDigraph, root::Int)::VVi

    # data
    V::Vi = collect(keys(intransitive_dag.V))
    A::Arcs = intransitive_dag.A

    # util
    adjList::Dict{Int, Vi} = Dict{Int, Vi}(map(i::Int -> i => Vi(), V))

    for (i::Int, j::Int) in A
        push!(adjList[i], j)
    end

    # output
    maximal_paths::VVi = VVi()

    # dfs
    nodes_buffer::Stack{Int} = Stack{Int}()

    function DFS(i::Int)
        
        # mark node
        push!(nodes_buffer, i)

        # edge case: save path
        isempty(adjList[i]) && push!(maximal_paths, Vi(reverse(collect(nodes_buffer))))

        # for neighbors
        for j::Int in adjList[i]
            DFS(j)
        end

        # unmark node
        pop!(nodes_buffer)
        
    end

    DFS(root)

    # return
    return maximal_paths
end

function filterMaximalPathsByWeight(intransitive_dag::InputDigraph, maximal_paths::VVi)::VVi
    return filter(
                  maximal_path::Vi -> sum(
                                          (i, j)::Tuple{Int, Int} -> intransitive_dag.distance[Arc(i, j)], 
                                          zip(maximal_path[1:end - 1], maximal_path[2:end])
                                         ) > 1, 
                  maximal_paths
                 )
end

function unpackPaths(maximal_paths::VVi, nodes_sccs::Dict{Int, Vi})::VVi
    return map(maximal_path::Vi -> reduce(union, map(i::Int -> nodes_sccs[i], maximal_path)), maximal_paths)
end

function getMaximalPathsAndInvalidArcs(data::SBRPData)::Tuple{VVi, ArcsSet}

    # get dominance digraph
    dominance_digraph::InputDigraph = getDominanceDigraph(data)

    # shrink digraph
    compact_dominance_digraph::InputDigraph, root::Int, nodes_sccs::Dict{Int, Vi} = getCompactDominanceDigraph(dominance_digraph)

    # transitive reduction
    intransitive_dag::InputDigraph = getDAGTransitiveReduction(compact_dominance_digraph, root)

    # enumerate maximal paths
    maximal_paths::VVi = getMaximalPathsDAG(intransitive_dag, root)

    # select those with at least two nodes
    filtered_maximal_paths::VVi = filterMaximalPathsByWeight(intransitive_dag, maximal_paths)

    # unpack paths
    maximal_unpacked_paths::VVi = unpackPaths(filtered_maximal_paths, nodes_sccs)

    # get preprocessed acs
    invalid_arcs::ArcsSet = ArcsSet()
    for (i::Int, j::Int) in dominance_digraph.A
        push!(invalid_arcs, Arc(i, j))
        push!(invalid_arcs, Arc(j, i))
    end

    # return
    return maximal_unpacked_paths, invalid_arcs
end

include("ip_model.jl")
include("cg_model.jl")
include("cp_model.jl")
include("concorde.jl")
