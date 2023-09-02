using MathOptInterface
using CPLEXCP
using ConstraintProgrammingExtensions

const MOI = MathOptInterface
const CP = ConstraintProgrammingExtensions

function runCompleteDigraphCPModel(data::SBRPData, app::Dict{String, Any})::Tuple{SBRPSolution, Dict{String, String}}

    # instance parameters
    
    @debug "Getting params"

    depot::Int                 = data.depot
    B::VVi                     = data.B
    A::Arcs                    = data.D.A
    T::Float64                 = data.T
    V::Vi                      = collect(keys(data.D.V))
    profits::Dict{Vi, Float64} = data.profits
    Vb::Si                     = getBlocksNodes(data)

    # aux data
    n::Int = length(V)

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
    t::Dict{Int, MOI.VariableIndex} = Dict{Int, MOI.VariableIndex}(map(i::Int -> i => first(MOI.add_constrained_variable(model, MOI.Interval(0.0, T))), V))

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
#    MOI.set(model, MOI.TimeLimitSec(), 3600)

    # Constraints
    
    @debug "Setup constraints"

    # No cycles
    
    @debug "No cycle constraints"

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
                       CP.DifferentFrom(depot),
                      )

    # Node visited
    
    @debug "Node visited"

    for i in V
        _ = MOI.add_constraint(
                           model,
                           MOI.VectorOfVariables([w[i], x[i]]),
                           CP.EquivalenceNot(MOI.EqualTo(1), MOI.EqualTo(i)),
                          )
    end

    # Block serviced
    
    @debug "Block serviced"

    for block::Vi in B

        visited_nodes::MOI.ScalarAffineFunction{Int} = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm{Int}.(
                                                                                                        vcat(map(i::Int -> 1, block), [-1]), 
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
        for j in V

            # edge case
            i == j && continue

            t_affine::MOI.ScalarAffineFunction{Float64} = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm{Float64}.([1, -1], [t[j], t[i]]), 
                                                                                   - time(data, Arc(i, j)))

            println(typeof(MOI.Utilities.operate(vcat, MOI.ScalarAffineFunction{Float64}, MOI.SingleVariable(x[i]), t_affine)))

            # if x[i] == j then t[j] = t[i] + time(i, j)
            _ = MOI.add_constraint(
                               model,
                               MOI.Utilities.operate(vcat, MOI.ScalarAffineFunction{Float64}, MOI.SingleVariable(x[i]), t_affine),
                               CP.Implication(MOI.EqualTo(j), MOI.EqualTo(0)),
                              )
        end
    end
    
    # Tour time
    
    @debug "Tour time"

    # Solve the model.
    
    @debug "Solve model"

    MOI.optimize!(model)
    
    # Check if the solution is optimum.
    @show MOI.get(model, MOI.TerminationStatus())
    
    # Get the solution
   
    @debug "Get solution"

    info::Dict{String, String} = Dict{String, String}()

    serviced_blocks::VVi = filter(block::Vi -> MOI.get(model, MOI.VariablePrimal(), y[block]) > .5, B)
    tour::Vi = Vi()

    # DFS
    tourAdjList::Dict{Int, Int} = Dict{Int, Int}(map(i::Int -> i => V[Int(MOI.get(model, MOI.VariablePrimal(), x[i]))], V))

    curr::Int = depot 
    push!(tour, curr)

    while true
        next::Int = tourAdjList[curr]
        push!(tour, next)  

        next == depot && break

        curr = next
    end

    solution::SBRPSolution = SBRPSolution(tour, serviced_blocks)

    println(tour)
#    for (k, v) in tourAdjList
#        println(k, " => ", v)
#    end
    for i in V
        if MOI.get(model, MOI.VariablePrimal(), w[i]) > .5
            println(i)
        end

    end

    # Return
   
    @debug "Return"

    return solution, info
end
