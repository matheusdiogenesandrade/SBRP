using MathOptInterface
using CPLEXCP
using ConstraintProgrammingExtensions

const MOI = MathOptInterface
const CP = ConstraintProgrammingExtensions

using Base.Threads

# add type 2 (CP focused)
function getIntersectionCutsCPFocused(data::SBRPData)::Arcs
    @debug "Obtaining intersection cuts"

    # data
    B::VVi = data.B 
    A::Arcs = data.D.A
    V::Vi = collect(keys(data.D.V))

    # output
    cuts2::ArcsSet = ArcsSet() # an arc (i, j) belongs to cuts2, whether i is visited than j should not be visited
    cliques::Set{VVi}   = Set{VVi}()

    # get nodes of each block
    Vb::Si = getBlocksNodes(data)
    nodes_blocks::Dict{Int, VVi} = Dict{Int, VVi}(map(i::Int -> i => filter(block::Vi -> i in block, B), Vi(collect(Vb))))

    # get cliques
    for i::Int in Vb

        @debug @sprintf("Obtaining cliques for the node %d", i)

        # data
        blocks::VVi     = nodes_blocks[i]
        blocks_num::Int = length(blocks)

        # edge case
        blocks_num == 1 && continue

        push!(cliques, blocks)
    end

    # get cuts
    for clique::VVi in cliques
        
        # edge case
        length(clique) <= 1 && continue

        @debug "Clique length $(length(clique))"

        # get intersection
        intersection::Vi = ∩(clique...)

        intersection_arcs::Arcs = reduce(vcat, map(i::Int -> δ⁺(A, i), intersection))

        @debug "Clique intersection $intersection"

        ######## intersection cuts 2 ########
        
        covered_nodes::Vi = filter(i::Int -> ⊆(nodes_blocks[i], clique) && !in(i, intersection), collect(Vb))

        for a::Arc in χ(intersection, covered_nodes)
            push!(cuts2, a)
        end

    end
    
    return collect(cuts2)
end

#=
# add type 2 (CP focused)
function getIntersectionCutsCPFocused(data::SBRPData)::Arcs

    # data
    Vb::Si = getBlocksNodes(data)
    B::VVi                     = data.B
    A::Arcs                    = data.D.A
    nodes_blocks::Dict{Int, VVi} = Dict{Int, VVi}(map(i::Int -> i => filter(block::Vi -> i in block, B), Vi(collect(Vb))))

    # output
    pairs::Arcs = Arcs()
    cliques::Vector{VVi} = Vector{VVi}()

    # get cliques
    for i in Vb

        @debug @sprintf("Obtaining cliques for the node %d", i)

        # data
        blocks::VVi = nodes_blocks[i]
        blocks_num::Int = length(blocks)

        # edge case
        blocks_num == 1 && continue

        # store
        for clique::VVi in collect(Combinatorics.powerset(blocks, 2))
            push!(cliques, clique)
        end

    end

    #=
    # max clique model
    max_clique::Model = direct_model(CPLEX.Optimizer())

    set_silent(max_clique)

    @variable(max_clique, z[block in B], Bin)

    @objective(max_clique, Max, ∑(z[block] for block in B))

    @constraint(max_clique, [(block, block′) in [(B[i], B[j]) for i in 1:length(B) for j in i + 1:length(B) if isempty(∩(B[i], B[j]))]], 1 >= z[block] + z[block′])
    @constraint(max_clique, ∑(z[block] for block in B) >= 2)

    while true
        # solve
        optimize!(max_clique)

        # base case
        termination_status(max_clique) == MOI.INFEASIBLE && break 

        # get clique
        B′::VVi = filter(block::Vi -> value(z[block]) > 0.5, B)

        # check if it is a clique 
        if all(
               (block′, block)::Pair{Vi, Vi} -> !isempty(∩(block, block′)), 
               filter((block′, block)::Pair{Vi, Vi} -> block′ ≠ block, χ(B, B′))
              ) 
            throw(InvalidStateException("It is not a clique"))
        end

        # store clique
        push!(cliques, B′)

        # update clique model
        @constraint(max_clique, sum(block::Vi -> z[block], B′) <= length(B′) - 1)
    end
    =#

    # get cuts
    for clique::VVi in cliques

        # edge case
        length(clique) <= 1 && continue

        @debug "Clique length $(length(clique))"

        # get intersection
        intersection::Vi = ∩(clique...)

        # get exclusive clique nodes
        exclusive_nodes::Vi = setdiff(∪(clique...), ∪(setdiff(B, clique)...))

        for (i::Int, j::Int) in χ(intersection, exclusive_nodes)
            i == j && continue
            push!(pairs, Pair{Int, Int}(i, j))
        end

    end

    return pairs
end
=#

function runCOPCompleteDigraphCPModel(data::SBRPData, app::Dict{String, Any})::Tuple{SBRPSolution, Dict{String, String}}

    # instance parameters
    
    @debug "Getting params"

    depot::Int                 = data.depot
    B::VVi                     = data.B
    A::Arcs                    = data.D.A
    V::Vi                      = sort(collect(keys(data.D.V))) # sorting for better debugging
    profits::Dict{Vi, Float64} = data.profits

    # aux data
    EPS = 10000;
    Vb::Si                      = getBlocksNodes(data)
    times::Dict{Arc, Int32}     = Dict{Arc, Int32}(map(a::Arc -> a => Int32(floor(time(data, a) * EPS)), A))
    T::Int32                    = Int32(floor(data.T * EPS))
    blockTimes::Dict{Vi, Int32} = Dict{Vi, Int32}(map(block::Vi -> block => Int32(floor(blockTime(data, block) * EPS)), B))

    nodes_blocks::Dict{Int, VVi} = Dict{Int, VVi}(map(i::Int -> i => filter(block::Vi -> i in block, B), Vi(collect(Vb))))

    n::Int = length(V)
    node_idx_map::Dict{Int, Int} = Dict{Int, Int}(map((idx, i)::Tuple{Int, Int} -> i => idx, enumerate(V)))

    # create model
    
    @debug "Create CP model"

    model = CPLEXCP.Optimizer()

    # variables
    
    @debug "Setup variables"

    # next node's index
    x::Dict{Int, MOI.VariableIndex} = Dict{Int, MOI.VariableIndex}(map(i::Int -> i => first(MOI.add_constrained_variable(model, MOI.Integer())), V))

    # node visited
    w::Dict{Int, MOI.VariableIndex} = Dict{Int, MOI.VariableIndex}(map(i::Int -> i => first(MOI.add_constrained_variable(model, MOI.ZeroOne())), V))

    # block flag
    y::Dict{Vi, MOI.VariableIndex} = Dict{Vi, MOI.VariableIndex}(map(block::Vi -> block => first(MOI.add_constrained_variable(model, MOI.ZeroOne())), B))

    # incurred time
    t::Dict{Int, MOI.VariableIndex} = Dict{Int, MOI.VariableIndex}(map(i::Int -> i => first(MOI.add_constrained_variable(model, MOI.Interval{Int32}(0, T))), V))

    # Domains
    
    @debug "Setup domains"

    idxs_set::Si = Si(1:n)

    # Nodes
    for i in V
        MOI.add_constraint(model, x[i], CP.Domain{Int}(idxs_set))
    end

    # Objective function
    
    @debug "Setup objective function"

    objective_function::MOI.ScalarAffineFunction{Int} = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(
                                                             map(block::Vi -> Int(profits[block]), B), 
                                                             map(block::Vi -> y[block], B)
                                                            ), 0)

    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    MOI.set(model, MOI.ObjectiveFunction{MathOptInterface.ScalarAffineFunction{Int}}(), objective_function)

    # Constraints
    
    @debug "Setup constraints"

    # Continuity constriants
    
    @debug "Continuity constraints"

    _ = MOI.add_constraint(
                       model,
                       MOI.VectorOfVariables(map(i::Int -> x[i], V)),
                       CP.AllDifferent(n),
                      )
    # Depot
    
    @debug "Depot constraint"

    _ = MOI.add_constraint(
                       model,
                       MOI.SingleVariable(x[depot]),
                       CP.DifferentFrom(node_idx_map[depot]),
                      )

#    _ = MOI.add_constraint(
#                           model,
#                           MOI.SingleVariable(t[depot]),
#                           MOI.EqualTo{Int32}(0)
#                          )
    CPLEXCP._set_ub(model, MOI.SingleVariable(t[depot]), MOI.EqualTo{Int32}(0).value)

    # Node visited
    
    @debug "Node visited"

    for (idx::Int, i::Int) in enumerate(V)
        _ = MOI.add_constraint(
                           model,
                           MOI.VectorOfVariables([w[i], x[i]]),
                           CP.EquivalenceNot(MOI.EqualTo(1), MOI.EqualTo(idx)),
                          )
    end

    # Block serviced
    
    @debug "Block serviced"

    for block::Vi in B

        visited_nodes::MOI.ScalarAffineFunction{Int} = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm{Int}.(
                                                                                                        vcat(ones(Int, length(block)),   [-1]), 
                                                                                                        vcat(map(i::Int -> w[i], block), [y[block]])
                                                                                                       ), 0)
        _ = MOI.add_constraint(
                               model,
                               visited_nodes,
                               MOI.GreaterThan(0)
                              )
    end

    # Incurred time

    @debug "Incurred time"

    for i in V
        # time 0
        # if x[i] == i then t[i] = 0.0
        implication::CP.Implication = CP.Implication(MOI.EqualTo{Int32}(node_idx_map[i]), MOI.EqualTo{Int32}(0))

        constr = CPLEXCP.cpo_java_imply(
                                        model.inner,
                                        CPLEXCP._build_constraint(model, x[i], implication.antecedent),
                                        CPLEXCP._build_constraint(model, t[i], implication.consequent),
                                       )
        CPLEXCP.cpo_java_add(model.inner, constr)

        # time limit
        # t[i] + \sum_{b \in B} y[b] * t_b  <= T
        affine = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm{Int32}.(
                                                                     vcat(map(block::Vi -> blockTimes[block], B), [1]), 
                                                                     vcat(map(block::Vi -> y[block], B), [t[i]])
                                                                    ), Int32(0))

        _ = MOI.add_constraint(
                               model,
                               affine,
                               MOI.LessThan(T)
                              )
       
    end

    # Tour time

    @debug "Tour time"

    for (i::Int, j::Int) in A

        # edge case
        j == depot && continue

        affine::MOI.ScalarAffineFunction{Int32} = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm{Int32}.([1, -1], [t[j], t[i]]), Int32(0))

        # if x[i] == j then t[j] = t[i] + time(i, j)
        implication::CP.Implication = CP.Implication(MOI.EqualTo{Int32}(node_idx_map[j]), MOI.EqualTo{Int32}(times[Arc(i, j)]))

        constr = CPLEXCP.cpo_java_imply(
                                        model.inner,
                                        CPLEXCP._build_constraint(model, x[i], implication.antecedent),
                                        CPLEXCP._build_constraint(model, affine, implication.consequent),
                                       )
        CPLEXCP.cpo_java_add(model.inner, constr)

    end

    # improvements

    @debug "Improvements"

    for block::Vi in B
        selected_nodes::Vi = filter(i::Int -> length(nodes_blocks[i]) == 1, block)
 
        visited_nodes::MOI.ScalarAffineFunction{Int} = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm{Int}.(
                                                                                                           vcat([1], map(_ -> -1, selected_nodes)),
                                                                                                           vcat([y[block]], map(i::Int -> w[i], selected_nodes)) 
                                                                                                          ), 0)
        _ = MOI.add_constraint(
                               model,
                               visited_nodes,
                               MOI.GreaterThan(0)
                              )
    end

    for block::Vi in B

        for (i::Int, j::Int) in χ(block, block)
            i == j && continue

            if length(nodes_blocks[i]) == 1 || length(nodes_blocks[j]) == 1
                _ = MOI.add_constraint(
                                       model,
                                       MOI.SingleVariable(x[i]),
                                       CP.DifferentFrom(node_idx_map[j]),
                                      )
            end

        end
    end

    # Intersection cuts

    @debug "Intersection cuts"

    # get intersection cuts
    intersection_cuts1::Vector{Arcs}, intersection_cuts2::Vector{Arcs} = getIntersectionCuts(data) 

    # add type 1
    #=
    for arcs::Arcs in intersection_cuts1 
        for (i::Int, j::Int) in arcs
            _ = MOI.add_constraint(
                                   model,
                                   MOI.SingleVariable(x[i]),
                                   CP.DifferentFrom(node_idx_map[j]),
                                  )
        end
    end
    =#

    # add type 2
#    for arcs::Arcs in intersection_cuts2
#        selected_nodes::Vi = collect(Si(map((i, j)::Pair{Int, Int} -> i, arcs)))
# 
#        visited_nodes::MOI.ScalarAffineFunction{Int} = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm{Int}.(
#                                                                                                           ones(Int, length(selected_nodes)), 
#                                                                                                           map(i::Int -> w[i], selected_nodes)
#                                                                                                          ), 0)
#        _ = MOI.add_constraint(
#                               model,
#                               visited_nodes,
#                               MOI.LessThan(1)
#                              )
#    end


    # add type 3
    pairs::Arcs = getIntersectionCutsCPFocused(data)

    #=
    for (i::Int, j::Int) in pairs
        implication::CP.Implication = CP.Implication(MOI.EqualTo{Int32}(1), MOI.EqualTo{Int32}(0))

        constr = CPLEXCP.cpo_java_imply(
                                        model.inner,
                                        CPLEXCP._build_constraint(model, w[i], implication.antecedent),
                                        CPLEXCP._build_constraint(model, w[j], implication.consequent),
                                       )
        CPLEXCP.cpo_java_add(model.inner, constr)
    end
    =#

    # parameters
    
    @debug "Setting parameters"
    
    cpo_java_setdoubleparameter(model.inner, "TimeLimit", 3600)
    cpo_java_setintparameter(model.inner, "Workers", Int32(nthreads()))

    # Solve the model.
    
    @debug "Solve model"
    elapsed_time::Float64 = @elapsed MOI.optimize!(model)
    
    # Check if the solution is optimum.
    @show MOI.get(model, MOI.TerminationStatus())
    
    # Get the solution
   
    @debug "Get solution"

    serviced_blocks::VVi = filter(block::Vi -> MOI.get(model, MOI.VariablePrimal(), y[block]) > .5, B)
    tour::Vi = Vi()

    # DFS
    tourAdjList::Dict{Int, Int} = Dict{Int, Int}(
                                                 map(i::Int -> i => V[Int(MOI.get(model, MOI.VariablePrimal(), x[i]))], V)
                                                )
    curr::Int = depot 
    push!(tour, curr)

    while true
        next::Int = tourAdjList[curr]
        push!(tour, next)  

        next == depot && break

        curr = next
    end

    solution::SBRPSolution = SBRPSolution(tour, serviced_blocks)

    cpo_java_release(model.inner)

    # Get the info
   
    @debug "Get info"

    info::Dict{String, String} = Dict{String, String}()
    info["cost"] = string(sum(block::Vi -> profits[block], serviced_blocks))
    info["solverTime"] = string(elapsed_time)


    # Return
   
    @debug "Return"

    return solution, info
end

function runCSPCompleteDigraphCPModel(data::SBRPData, app::Dict{String, Any})::Tuple{SBRPSolution, Dict{String, String}}

    # instance parameters
    
    @debug "Getting params"

    depot::Int                 = data.depot
    B::VVi                     = data.B
    A::Arcs                    = data.D.A
    V::Vi                      = sort(collect(keys(data.D.V))) # sorting for better debugging
    profits::Dict{Vi, Float64} = data.profits

    # aux data
    EPS::Int                    = all(value::Float64 -> value == floor(value), map(a::Arc -> time(data, a), A)) ? 1 : 10000
    Vb::Si                      = getBlocksNodes(data)
    times::Dict{Arc, Int32}     = Dict{Arc, Int32}(map(a::Arc -> a => Int32(floor(time(data, a) * EPS)), A))
    T::Int32                    = Int32(floor(data.T * EPS))

    nodes_blocks::Dict{Int, VVi} = Dict{Int, VVi}(map(i::Int -> i => filter(block::Vi -> i in block, B), Vi(collect(Vb))))

    n::Int = length(V)
    node_idx_map::Dict{Int, Int} = Dict{Int, Int}(map((idx, i)::Tuple{Int, Int} -> i => idx, enumerate(V)))

    # create model
    
    @debug "Create CP model"

    model = CPLEXCP.Optimizer()

    # variables
    
    @debug "Setup variables"

    # next node's index
    x::Dict{Int, MOI.VariableIndex} = Dict{Int, MOI.VariableIndex}(map(i::Int -> i => first(MOI.add_constrained_variable(model, MOI.Integer())), V))

    # node visited
    w::Dict{Int, MOI.VariableIndex} = Dict{Int, MOI.VariableIndex}(map(i::Int -> i => first(MOI.add_constrained_variable(model, MOI.ZeroOne())), V))

    # incurred time
    t::Dict{Int, MOI.VariableIndex} = Dict{Int, MOI.VariableIndex}(map(i::Int -> i => first(MOI.add_constrained_variable(model, MOI.Interval{Int32}(0, T))), V))

    # Domains
    
    @debug "Setup domains"

    idxs_set::Si = Si(1:n)

    # Nodes
    for i::Int in V
        MOI.add_constraint(model, x[i], CP.Domain{Int}(idxs_set))
    end

    # Objective function
    
    @debug "Setup objective function"

    objective_function::MOI.ScalarAffineFunction{Int} = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(
                                                             map(i::Int -> 1, V), 
                                                             map(i::Int -> times[Arc(i, x[i])], V)
                                                            ), 0)

    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    MOI.set(model, MOI.ObjectiveFunction{MathOptInterface.ScalarAffineFunction{Int}}(), objective_function)

    # Constraints
    
    @debug "Setup constraints"

    # Continuity constriants
    
    @debug "Continuity constraints"

    _ = MOI.add_constraint(
                       model,
                       MOI.VectorOfVariables(map(i::Int -> x[i], V)),
                       CP.AllDifferent(n),
                      )
    # Depot
    
    @debug "Depot constraint"

    _ = MOI.add_constraint(
                       model,
                       MOI.SingleVariable(x[depot]),
                       CP.DifferentFrom(node_idx_map[depot]),
                      )

#    _ = MOI.add_constraint(
#                           model,
#                           MOI.SingleVariable(t[depot]),
#                           MOI.EqualTo{Int32}(0)
#                          )
    CPLEXCP._set_ub(model, MOI.SingleVariable(t[depot]), MOI.EqualTo{Int32}(0).value)

    # Node visited
    
    @debug "Node visited"

    for (idx::Int, i::Int) in enumerate(V)
        _ = MOI.add_constraint(
                           model,
                           MOI.VectorOfVariables([w[i], x[i]]),
                           CP.EquivalenceNot(MOI.EqualTo(1), MOI.EqualTo(idx)),
                          )
    end

    # Block serviced
    
    @debug "Block serviced"

    for block::Vi in B

        visited_nodes::MOI.ScalarAffineFunction{Int} = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm{Int}.(
                                                                                                        vcat(ones(Int, length(block)),   [-1]), 
                                                                                                        vcat(map(i::Int -> w[i], block), [y[block]])
                                                                                                       ), 0)
        _ = MOI.add_constraint(
                               model,
                               visited_nodes,
                               MOI.GreaterThan(0)
                              )
    end

    # Incurred time

    @debug "Incurred time"

    for i in V
        # time 0
        # if x[i] == i then t[i] = 0.0
        implication::CP.Implication = CP.Implication(MOI.EqualTo{Int32}(node_idx_map[i]), MOI.EqualTo{Int32}(0))

        constr = CPLEXCP.cpo_java_imply(
                                        model.inner,
                                        CPLEXCP._build_constraint(model, x[i], implication.antecedent),
                                        CPLEXCP._build_constraint(model, t[i], implication.consequent),
                                       )
        CPLEXCP.cpo_java_add(model.inner, constr)

        # time limit
        # t[i] + \sum_{b \in B} y[b] * t_b  <= T
        affine = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm{Int32}.(
                                                                     vcat(map(block::Vi -> blockTimes[block], B), [1]), 
                                                                     vcat(map(block::Vi -> y[block], B), [t[i]])
                                                                    ), Int32(0))

        _ = MOI.add_constraint(
                               model,
                               affine,
                               MOI.LessThan(T)
                              )
       
    end

    # Tour time

    @debug "Tour time"

    for (i::Int, j::Int) in A

        # edge case
        j == depot && continue

        affine::MOI.ScalarAffineFunction{Int32} = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm{Int32}.([1, -1], [t[j], t[i]]), Int32(0))

        # if x[i] == j then t[j] = t[i] + time(i, j)
        implication::CP.Implication = CP.Implication(MOI.EqualTo{Int32}(node_idx_map[j]), MOI.EqualTo{Int32}(times[Arc(i, j)]))

        constr = CPLEXCP.cpo_java_imply(
                                        model.inner,
                                        CPLEXCP._build_constraint(model, x[i], implication.antecedent),
                                        CPLEXCP._build_constraint(model, affine, implication.consequent),
                                       )
        CPLEXCP.cpo_java_add(model.inner, constr)

    end

    # improvements

    @debug "Improvements"

    for block::Vi in B
        selected_nodes::Vi = filter(i::Int -> length(nodes_blocks[i]) == 1, block)
 
        visited_nodes::MOI.ScalarAffineFunction{Int} = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm{Int}.(
                                                                                                           vcat([1], map(_ -> -1, selected_nodes)),
                                                                                                           vcat([y[block]], map(i::Int -> w[i], selected_nodes)) 
                                                                                                          ), 0)
        _ = MOI.add_constraint(
                               model,
                               visited_nodes,
                               MOI.GreaterThan(0)
                              )
    end

    for block::Vi in B

        for (i::Int, j::Int) in χ(block, block)
            i == j && continue

            if length(nodes_blocks[i]) == 1 || length(nodes_blocks[j]) == 1
                _ = MOI.add_constraint(
                                       model,
                                       MOI.SingleVariable(x[i]),
                                       CP.DifferentFrom(node_idx_map[j]),
                                      )
            end

        end
    end

    # Intersection cuts

    @debug "Intersection cuts"

    # get intersection cuts
    intersection_cuts1::Vector{Arcs}, intersection_cuts2::Vector{Arcs} = getIntersectionCuts(data) 

    # add type 1
    #=
    for arcs::Arcs in intersection_cuts1 
        for (i::Int, j::Int) in arcs
            _ = MOI.add_constraint(
                                   model,
                                   MOI.SingleVariable(x[i]),
                                   CP.DifferentFrom(node_idx_map[j]),
                                  )
        end
    end
    =#

    # add type 2
#    for arcs::Arcs in intersection_cuts2
#        selected_nodes::Vi = collect(Si(map((i, j)::Pair{Int, Int} -> i, arcs)))
# 
#        visited_nodes::MOI.ScalarAffineFunction{Int} = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm{Int}.(
#                                                                                                           ones(Int, length(selected_nodes)), 
#                                                                                                           map(i::Int -> w[i], selected_nodes)
#                                                                                                          ), 0)
#        _ = MOI.add_constraint(
#                               model,
#                               visited_nodes,
#                               MOI.LessThan(1)
#                              )
#    end


    # add type 3
    pairs::Arcs = getIntersectionCutsCPFocused(data)

    #=
    for (i::Int, j::Int) in pairs
        implication::CP.Implication = CP.Implication(MOI.EqualTo{Int32}(1), MOI.EqualTo{Int32}(0))

        constr = CPLEXCP.cpo_java_imply(
                                        model.inner,
                                        CPLEXCP._build_constraint(model, w[i], implication.antecedent),
                                        CPLEXCP._build_constraint(model, w[j], implication.consequent),
                                       )
        CPLEXCP.cpo_java_add(model.inner, constr)
    end
    =#

    # parameters
    
    @debug "Setting parameters"
    
    cpo_java_setdoubleparameter(model.inner, "TimeLimit", 3600)
    cpo_java_setintparameter(model.inner, "Workers", Int32(nthreads()))

    # Solve the model.
    
    @debug "Solve model"
    elapsed_time::Float64 = @elapsed MOI.optimize!(model)
    
    # Check if the solution is optimum.
    @show MOI.get(model, MOI.TerminationStatus())
    
    # Get the solution
   
    @debug "Get solution"

    serviced_blocks::VVi = filter(block::Vi -> MOI.get(model, MOI.VariablePrimal(), y[block]) > .5, B)
    tour::Vi = Vi()

    # DFS
    tourAdjList::Dict{Int, Int} = Dict{Int, Int}(
                                                 map(i::Int -> i => V[Int(MOI.get(model, MOI.VariablePrimal(), x[i]))], V)
                                                )
    curr::Int = depot 
    push!(tour, curr)

    while true
        next::Int = tourAdjList[curr]
        push!(tour, next)  

        next == depot && break

        curr = next
    end

    solution::SBRPSolution = SBRPSolution(tour, serviced_blocks)

    cpo_java_release(model.inner)

    # Get the info
   
    @debug "Get info"

    info::Dict{String, String} = Dict{String, String}()
    info["cost"] = string(sum(block::Vi -> profits[block], serviced_blocks))
    info["solverTime"] = string(elapsed_time)


    # Return
   
    @debug "Return"

    return solution, info
end
