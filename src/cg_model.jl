using Printf

mutable struct BBNode
   id::Int
   model::Model
   pred::Union{BBNode, Nothing}
   children::Vector{BBNode}
   BBNode(id::Int, model::Model, pred::Union{BBNode, Nothing}) = new(id, deepcopy(model), pred)
end

#=
Get subtour cuts
input: 
- V::Vi
- A::Arcs
- model::Model is a Mathematical Programming model
output:
- components::Vector{Tuple{Vi, Int, Int}} is the set of sets of nodes to separate.
=# 
function getSubtourCuts(V::Vi, A::Arcs, model::Model, app::Dict{String, Any})::Vector{Tuple{Vi, Int, Int}}

    @debug "Obtaining subtour cuts"

    # setup
    x = model[:x]

    # outputs
    nodes_sets::Set{Tuple{Si, Int, Int}} = Set{Tuple{Si, Int, Int}}()

    # helpers
    V_map::Dict{Int, Int} = Dict{Int, Int}(
                        map(
                            (idx, i)::Tuple{Int, Int} -> i => idx, 
                            enumerate(V)
                           )
                       )
    V_map_reverse::Dict{Int, Int} = Dict{Int, Int}(
                                                   map(
                                                       (idx, i)::Tuple{Int, Int} -> idx => i, 
                                                       enumerate(V)
                                                      )
                                                  )

    n::Int = length(V)

    # get values
    x_val::ArcCostMap = ArcCostMap(map(a::Arc -> a => value(x[a]), A))
    z_val::Dict{Int, Float64} = Dict{Int, Float64}(map(i::Int -> i => sum(a::Arc -> x_val[a], δ⁺(A, i)), V))

    # get subsets
    g = SparseMaxFlowMinCut.ArcFlow[]
    M::Int = 100000

    # used arcs
    A′::Arcs = filter(a::Arc -> x_val[a] > EPS, A)

    # used nodes
    V′::Vi = filter(i::Int -> z_val[i] > EPS, V)

    # mounting graph
    for (i::Int, j::Int) in A′
        push!(g, SparseMaxFlowMinCut.ArcFlow(V_map[i], V_map[j], trunc(floor(x_val[Arc(i, j)], digits=5) * M)))
    end

    # get subsets
    for source::Int in V′
        for target::Int in V′

            # edge case 
            source == target && continue

            # init
            maxFlow::Float64, flows, set = SparseMaxFlowMinCut.find_maxflow_mincut(
                                                                                   SparseMaxFlowMinCut.Graph(n, g), 
                                                                                   V_map[source], 
                                                                                   V_map[target]
                                                                                  )
            flow::Float64 = maxFlow / M

            # base case: In the same component
            set[V_map[target]] == 1 && continue

            # base case: Condition not met
            flow + 1e-2 >= z_val[source] + z_val[target] - 1 && continue

            # get set
            S::Si = Si(
                       map(
                           i::Int -> V_map_reverse[i], 
                           filter(i::Int -> set[i] == 1, 1:n)
                          )
                      )

            # store
            entry::Tuple{Si, Int, Int} = Tuple{Si, Int, Int}((S, source, target))
            push!(nodes_sets, entry)
            
        end
    end

    return Vector{Tuple{Vi, Int, Int}}(
                                       map(
                                           (S, i, j)::Tuple{Si, Int, Int} -> Tuple{Vi, Int, Int}((collect(S), i, j)), 
                                           collect(nodes_sets)
                                          )
                                      )
end

#=
Build SBRP model
input: 
- data::SBRPData is a SBRP instance
- app::Dict{String, Any} is the relation of application parameters
=# 
function runCOPColumnGenerationModel(data::SBRPData, app::Dict{String, Any})::Tuple{SBRPSolution, Dict{String, String}}

    # instance parameters
    depot::Int                 = data.depot
    B::VVi                     = data.B
    A::Arcs                    = data.D.A
    T::Float64                 = data.T
    V::Dict{Int, Vertex}       = data.D.V
    profits::Dict{Vi, Float64} = data.profits
    Vb::Si                     = getBlocksNodes(data)
    nodes_blocks::Dict{Int, VVi} = Dict{Int, VVi}(map(i::Int -> i => filter(block::Vi -> i in block, B), collect(Vb)))

    # nodes indexes
    V_seq::Vi = collect(keys(V))
    V_idx::Dict{Int, Int} = Dict{Int, Int}(i => idx for (idx::Int, i::Int) in enumerate(V_seq))

    # output
    info::Dict{String, String} = Dict{String, String}("lazyCuts" => "0") 

    # MTZ data
#    A′::Arcs = filter((i, j)::Arc -> j != depot, A)
#    V′::Si = setdiff(keys(V), depot)
#    P′::Arcs = filter((i, j)::Arc -> i < j, χ(V′))

    # Initial arcs
    m::Int = length(B)
    A′′::Arcs = []
    for block::Vi in B
        push!(A′′, Arc(depot, first(block)))
        push!(A′′, Arc(first(block), depot))
    end
    for (block1::Vi, block2::Vi) in zip(B[1:m - 1], B[2:m])
        push!(A′′, Arc(first(block1), first(block2)))
    end

    V′′::Vi = collect(filter(i::Int -> any(a::Arc -> a[1] == i || a[2] == i, A′′), keys(V)))

    # helpers

    function getDual(a::Arc, model::Model)::Float64
        # params
        i::Int, j::Int = a
        selected_blocks::VVi = i == depot ? [] : nodes_blocks[i]

        i_key::JuMP.Containers.DenseAxisArrayKey{Tuple{Int64}} = JuMP.Containers.DenseAxisArrayKey{Tuple{Int64}}((i,))
        j_key::JuMP.Containers.DenseAxisArrayKey{Tuple{Int64}} = JuMP.Containers.DenseAxisArrayKey{Tuple{Int64}}((j,))

        degree_i_dual::Float64              = i_key in keys(model[:degree]) ? - dual(model[:degree][i_key]) : 0.0
        degree_j_dual::Float64              = j_key in keys(model[:degree]) ? dual(model[:degree][j_key]) : 0.0
        degree_at_most_once_i_dual::Float64 = i_key in keys(model[:degree_at_most_once]) && i != depot ? - dual(model[:degree_at_most_once][i_key]) : 0.0
        block_lift_dual::Float64            = sum(block::Vi -> dual(model[:block_lift][block]), selected_blocks; init = 0.0)
        time_limit_dual::Float64            = - dual(model[:time_limit]) * time(data, a)
        mtz_ub2_i_dual::Float64             = i_key in keys(model[:mtz_ub2]) && i != depot ? dual(model[:mtz_ub2][i_key]) * T : 0.0
        mtz_ub3_j_dual::Float64             = j_key in keys(model[:mtz_ub3]) && j != depot ? dual(model[:mtz_ub3][j_key]) * T : 0.0

        #
        reduced_cost::Float64 = degree_i_dual + degree_j_dual + degree_at_most_once_i_dual + block_lift_dual + time_limit_dual + mtz_ub2_i_dual + mtz_ub3_j_dual
#        reduced_cost::Float64 = degree_i_dual + degree_j_dual + degree_at_most_once_i_dual + block_lift_dual + dual(model[:time_limit])

        if i == depot
            reduced_cost += dual(model[:depot_visit])
        end

        return - reduced_cost
    end

    function createModel(relax_x::Bool = false, relax_y::Bool = false, verbose::Bool = false)::Model

        # params
        A′::Arcs = filter((i, j)::Arc -> j != depot, A′′)
        V′::Vi = setdiff(V′′, depot)

        local model

        model::Model = direct_model(CPLEX.Optimizer())
        set_parameters(model, "CPX_PARAM_TILIM" => 3600)

        !verbose && set_silent(model)

        if relax_x
            @variable(model, x[a in A′′], lower_bound = 0, upper_bound = 1)
        else
            @variable(model, x[a in A′′], Bin)
        end
        if relax_y
            @variable(model, y[b in B], lower_bound = 0, upper_bound = 1)
        else
            @variable(model, y[b in B], Bin)
        end

        @objective(model, Max, sum(block::Vi -> data.profits[block] * y[block], B))

        @constraint(model, degree[i::Int in V′′], sum(a::Arc -> x[a], δ⁻(A′′, i); init = 0) == sum(a::Arc -> x[a], δ⁺(A′′, i); init = 0))
        @constraint(model, degree_at_most_once[i::Int in V′], sum(a::Arc -> x[a], δ⁺(A′′, i)) <= 1)
        @constraint(model, depot_visit, sum(a::Arc -> x[a], δ⁺(A′′, depot)) == 1)
        @constraint(model, block_lift[block::Vi in B], y[block] <= sum(a::Arc -> x[a], δ⁺(A′′, block)))
        @constraint(model, time_limit, sum(a::Arc -> time(data, a) * x[a], A′′) <= T - sum(block::Vi -> y[block] * blockTime(data, block), B))

        #=
        @constraint(model, degree[i::Int in V′′], sum(a::Arc -> x[a], δ⁻(A′′, i); init = 0) - sum(a::Arc -> x[a], δ⁺(A′′, i); init = 0) == 0)
        @constraint(model, degree_at_most_once[i::Int in V′′], sum(a::Arc -> x[a], δ⁺(A′′, i)) - 1 <= 0)
        @constraint(model, depot_visit, sum(a::Arc -> x[a], δ⁺(A′′, depot)) - 1 == 0)
        @constraint(model, block_lift[block::Vi in B], y[block] - sum(a::Arc -> x[a], δ⁺(A′′, block)) <= 0)
        @constraint(model, time_limit, sum(a::Arc -> time(data, a) * x[a], A′′) + sum(block::Vi -> y[block] * blockTime(data, block), B) - T <= 0)
        =#

        # improvements
#        @constraint(model, block3[block::Vi in B], y[block] - sum(x[a] for i in block for a in δ⁺(A′′, i) if length(nodes_blocks[i]) == 1) >= 0)
#        @constraint(model, block4[block::Vi in B], sum((i, j)::Arc -> x[(i, j)], filter((i, j)::Arc -> i in block && j in block && nodes_blocks[i] == nodes_blocks[j] && length(nodes_blocks[j]) == 1, A)) == 0)
#        @constraint(model, block4[block::Vi in B], sum(a::Arc -> x[a], filter((i, j)::Arc -> i in block && j in block && (length(nodes_blocks[i]) == 1 || length(nodes_blocks[j]) == 1), A′′), init = 0.0) == 0)
#        @constraint(model, subcycle_size_two[a::Arc in P′], x[a] + x[reverse(a)] <= 1)

        # MTZs
        if app["arcs-mtz"] # Arcs MTZ
        else # Nodes MTZ

#            @variable(model, t[i::Int in keys(V)], lower_bound = 0, upper_bound = T)
            @variable(model, t[i::Int in V′′], lower_bound = 0, upper_bound = T)

            @constraint(model, time_at_depot, t[depot] == 0.0)
#            @constraint(model, mtz[a::Arc in A′′], t[last(a)] >= t[first(a)] + x[a] * time(data, a) - (1 - x[a]) * T - x[reverse(a)] * time(data, reverse(a)))
            @constraint(model, mtz[a::Arc in A′],     t[last(a)] >= t[first(a)] + x[a] * time(data, a) - (1 - x[a]) * T)
            @constraint(model, mtz_ub1[i::Int in V′], t[i] <= T - sum(block::Vi -> y[block] * blockTime(data, block), B))
            @constraint(model, mtz_ub2[i::Int in V′], t[i] <= sum(a::Arc -> x[a], δ⁺(A′′, i)) * T)
            @constraint(model, mtz_ub3[i::Int in V′], t[i] <= sum(a::Arc -> x[a], δ⁻(A′′, i)) * T)
            #=
            @constraint(model, mtz[a::Arc in A′′],     t[first(a)] + x[a] * time(data, a) - (1 - x[a]) * T - t[last(a)] <= 0)
            @constraint(model, mtz_ub1[i::Int in V′′], t[i] - (T - sum(block::Vi -> y[block] * blockTime(data, block), B)) <= 0)
            @constraint(model, mtz_ub2[i::Int in V′′], t[i] - sum(a::Arc -> x[a], δ⁺(A′′, i)) * T <= 0)
            @constraint(model, mtz_ub3[i::Int in V′′], t[i] - sum(a::Arc -> x[a], δ⁻(A′′, i)) * T <= 0)
            =#
        end

        # lb
        if app["lb"] != nothing 
            lb::Int = parse(Int, app["lb"])
            @constraint(model, sum(block::Vi -> data.profits[block] * y[block], B) >= lb)
        end

        # warm start
        if app["warm-start-solution"] != nothing 

            # get route
            (route::Vi, B_selected::VVi) = readSolution(app["warm-start-solution"], data)

            # get visited arcs
            route_arcs::ArcsSet = ArcsSet(zip(route[begin:end - 1], route[2:end]))

            # get non visited arcs
            non_visited_arcs::Arcs = filter(a::Arc -> !in(a, route_arcs), A)

            # fix variables
            for a::Arc in route_arcs
                set_start_value(x[a], 1.0)
            end
            for a::Arc in non_visited_arcs
                set_start_value(x[a], 0.0)
            end

        end

        # registries
        # variables
        model[:x], model[:y] = x, y

        # constraints
        model[:depot_visit]         = depot_visit
        model[:degree]              = degree 
        model[:degree_at_most_once] = degree_at_most_once
        model[:time_limit]          = time_limit
        model[:block_lift]          = block_lift
        # mtz
        model[:time_at_depot] = time_at_depot
        model[:mtz_ub1] = mtz_ub1
        model[:mtz_ub2] = mtz_ub2
        model[:mtz_ub3] = mtz_ub3

        return model
    end

    k::Int = 0
    K::Int = 1000
    while k < K

        local model

        # creating model with both variables (x and y) relaxed
#        model::Model = createModel(true, true)
        model::Model = createModel(true, true)

        optimize!(model)

        # filter arcs
        selected_arcs::Arcs = filter(a::Arc -> !in(a, A′′), A)

        # edge case
        isempty(selected_arcs) && break

        # get best arc
#        a::Arc = argmin(a::Arc -> getDual(a, model), selected_arcs)
        a::Arc = first(selected_arcs)
        best_dual::Float64 = getDual(a, model)

        for a_::Arc in selected_arcs 
            curr_dual::Float64 = getDual(a_, model) 

            if curr_dual > best_dual
                a = a_
                best_dual = curr_dual
            end
        end

        # edge case
        if best_dual <= 0

            @debug @sprintf("Best dual (%f) is negative or null , then break", best_dual)

            break
        end

        # store arc
        i::Int, j::Int = a
        push!(A′′, a)
        !in(i, V′′) && push!(V′′, i)
        !in(j, V′′) && push!(V′′, j)
        (j != depot && !in(Arc(j, depot), A′′)) && push!(A′′, Arc(j, depot))
        (depot != i && !in(Arc(depot, i), A′′)) && push!(A′′, Arc(depot, i))

        # log
        println(@sprintf("Iteration %d", k))
        println(@sprintf("\tObjective function: %f", objective_value(model)))
        println(@sprintf("\tBest arc: (%d, %d)", V_idx[i], V_idx[j]))
        println(@sprintf("\tBest dual: %f", best_dual))

        # increment
        k += 1
    end

    # creating model with both variables (x and y) relaxed
    model::Model = createModel(false, false, true)

    optimize!(model)

    println("Objective function: ", objective_value(model))
    println("Relative GAP: ", relative_gap(model))

    for a::Arc in A′′
        i, j = a
        if value(model[:x][a]) > 0.5
            println(V_idx[i], ", ", V_idx[j])
        end
    end

    for block::Vi in B
        if value(model[:y][block]) > 0.5
            println(block)
        end
    end


    #=
    # getting initial relaxation with both variables <x and y> relaxed
   
    @debug "Getting initial relaxation"

    # creating model with both variables (x and y) relaxed
    model::Model = createModel(true, true)

    optimize!(model)

    info["initialLP"] = string(objective_value(model))

    # getting initial relaxation with only x relaxed (y integer)
    
    @debug "Getting initial relaxation with y as integer"

    model = createModel(true)

    info["yLPTime"] = string(@elapsed optimize!(model))

    info["yLP"] = string(objective_value(model))

    # integer model
#    model = createModel()
#    MOI.set(model, MOI.NumberOfThreads(), 1) # thread numbers

    # run
    info["solverTime"]  = string(@elapsed optimize!(model))
    info["cost"]        = @sprintf("%.2f", objective_value(model))
    info["relativeGAP"] = string(relative_gap(model))
    info["nodeCount"]   = string(node_count(model))

    # retrieve solution
    x, y = model[:x], model[:y]

    solution_arcs::Arcs = Arcs(filter(a::Arc -> value(x[a]) > 0.5, A))
    solution_blocks::VVi = VVi(filter(block::Vi -> value(y[block]) > 0.5, B))
    tour::Vi = Vi([depot])

    while !isempty(solution_arcs)

        arc_idx::Union{Int, Nothing} = findfirst(a::Arc -> first(a) == last(tour), solution_arcs)

        if arc_idx == nothing 
            error("Desired arc not found")
        end

        a::Arc = solution_arcs[arc_idx]

        i::Int, j::Int = first(a), last(a)

        push!(tour, j)

        deleteat!(solution_arcs, arc_idx)
    end

    solution::SBRPSolution   = SBRPSolution(tour, solution_blocks)

    # return
    return solution, info
    =#
    
end

#=
Build CSP model
input: 
- data::SBRPData is a SBRP instance
- app::Dict{String, Any} is the relation of application parameters
=# 
function runCSPColumnGenerationModel(
        data::SBRPData, 
        app::Dict{String, Any},
        time::Function = time
    )::Tuple{SBRPSolution, Dict{String, String}}

    # instance parameters
    depot::Int                   = data.depot
    B::VVi                       = data.B
    A::Arcs                      = data.D.A
    A_set::ArcsSet               = ArcsSet(data.D.A)
    T::Float64                   = data.T
    V::Vi                        = collect(keys(data.D.V))
    Vb::Si                       = getBlocksNodes(data)
    nodes_blocks::Dict{Int, VVi} = Dict{Int, VVi}(map(i::Int -> i => filter(block::Vi -> i in block, B), collect(Vb)))
    m::Int                       = length(B)

    # nodes indexes
    V_seq::Vi = collect(V)
    V_idx::Dict{Int, Int} = Dict{Int, Int}(i => idx for (idx::Int, i::Int) in enumerate(V_seq))

    # output
    info::Dict{String, String} = Dict{String, String}("lazyCuts" => "0") 

    # Initial arcs
    A′′::Arcs = map(a::Tuple{Int, Int} -> Arc(a[1], a[2]), zip(V[begin:end - 1], V[begin + 1:end]))
    push!(A′′, Arc(last(V), first(V)))

    for (block1::Vi, block2::Vi) in zip(B[1:m - 1], B[2:m])

        i::Int, j::Int = first(block1), first(block2)
        a::Arc = Arc(i, j)

        if i != j && !(a in A′′)
            push!(A′′, a)
        end

    end

    # Subtour components
    subtour_components::Vector{Tuple{Vi, Int, Int}} = Vector{Tuple{Vi, Int, Int}}()

    if first(last(B)) != first(first(B)) && !(Arc(first(last(B)), first(first(B))) in A′′)
        push!(A′′, Arc(first(last(B)), first(first(B))))
    end

    # create model
    model::Model = direct_model(CPLEX.Optimizer())

    # params
    set_silent(model)
    # set_parameters(model, "CPX_PARAM_TILIM" => 3600)
    set_parameters(model, "CPX_PARAM_PREIND" => 0)
    # set_parameters(model, "CPXPARAM_MIP_Strategy_Search" => 1)
    set_parameters(model, "CPXPARAM_MIP_Cuts_FlowCovers" => -1)
    set_parameters(model, "CPXPARAM_MIP_Cuts_MIRCut" => -1)
    set_parameters(model, "CPXPARAM_MIP_Cuts_ZeroHalfCut" => -1)
    set_parameters(model, "CPXPARAM_MIP_Cuts_Gomory" => -1)
    set_parameters(model, "CPXPARAM_MIP_Cuts_Implied" => -1)
    set_parameters(model, "CPXPARAM_MIP_Cuts_LiftProj" => -1)
    set_parameters(model, "CPXPARAM_MIP_Cuts_Cliques" => -1)
    set_parameters(model, "CPXPARAM_MIP_Cuts_Covers" => -1)
    set_parameters(model, "CPXPARAM_MIP_Cuts_Disjunctive" => -1)
    set_parameters(model, "CPXPARAM_MIP_Cuts_PathCut" => -1)
    set_parameters(model, "CPXPARAM_MIP_Cuts_GUBCovers" => -1)
    set_parameters(model, "CPXPARAM_MIP_Limits_CutPasses" => -1)
    #
    set_parameters(model, "CPXPARAM_Preprocessing_Presolve" => 0)
    set_parameters(model, "CPXPARAM_Preprocessing_Aggregator" => 0)
    set_parameters(model, "CPXPARAM_Preprocessing_Relax" => 0)
    set_parameters(model, "CPXPARAM_Preprocessing_NumPass" => 0)
    set_parameters(model, "CPXPARAM_Preprocessing_RepeatPresolve" => 0)
    set_parameters(model, "CPXPARAM_Preprocessing_BoundStrength" => 0)
    set_parameters(model, "CPXPARAM_Preprocessing_CoeffReduce" => 0)
    set_parameters(model, "CPXPARAM_Preprocessing_Reduce" => 0)
    set_parameters(model, "CPXPARAM_Preprocessing_Symmetry" => 0)
    set_parameters(model, "CPXPARAM_Preprocessing_Fill" => 0)
    set_parameters(model, "CPXPARAM_Preprocessing_Dual" => 0)
    
    # setup
    @variable(model, x[a in copy(A′′)], lower_bound = 0, upper_bound = 1)
    model[:x] = x

    @objective(model, Min, sum(a::Arc -> time(data, a) * x[a], A′′))

    @constraint(model, degree[i::Int in V], sum(a::Arc -> x[a], δ⁻(A′′, i)) - sum(a::Arc -> x[a], δ⁺(A′′, i)) == 0)
    @constraint(model, degree_at_most_once[i::Int in V], 1 - sum(a::Arc -> x[a], δ⁺(A′′, i)) >= 0)
    @constraint(model, block_out_degree[block::Vi in B], sum(a::Arc -> x[a], δ⁺(A′′, block)) - 1 >= 0)
    @constraint(model, subtour[(S, i, j)::Tuple{Vi, Int, Int} in copy(subtour_components)], sum(a::Arc -> x[a], δ⁺(A′′, S)) - (sum(a::Arc -> x[a], δ⁺(A′′, i)) + sum(a::Arc -> x[a], δ⁺(A′′, j)) - 1) >= 0)

    # util
    function getDual(a::Arc)::Float64

        # params
        i::Int, j::Int = a
        selected_blocks::VVi = nodes_blocks[i]

        i_key::JuMP.Containers.DenseAxisArrayKey{Tuple{Int64}} = JuMP.Containers.DenseAxisArrayKey{Tuple{Int64}}((i,))
        j_key::JuMP.Containers.DenseAxisArrayKey{Tuple{Int64}} = JuMP.Containers.DenseAxisArrayKey{Tuple{Int64}}((j,))

        degree_i_dual::Float64              = - dual(degree[i_key])
        degree_j_dual::Float64              = dual(degree[j_key])
        degree_at_most_once_i_dual::Float64 = - dual(degree_at_most_once[i_key])
        block_out_degree_dual::Float64      = sum(block::Vi -> dual(block_out_degree[block]), filter(block::Vi -> a in δ⁺(A, block), selected_blocks); init = 0.0)
        subtour_dual::Float64               = 0.0

        for entry::Tuple{Vi, Int, Int} in subtour_components 
            S::Vi, i_::Int, j_::Int = entry
            if a in δ⁺(A′′, S)
                subtour_dual += dual(subtour[entry]) 
            end

            if i == i_ || i == j_
                subtour_dual -= dual(subtour[entry]) 
            end
        end

        #
        reduced_cost::Float64 = time(data, a) - (degree_i_dual + degree_j_dual + degree_at_most_once_i_dual + block_out_degree_dual + subtour_dual)

        return reduced_cost
    end

    function appendVariable(a::Arc)

        # store arc
        i::Int, j::Int = a

        # Create a new column
        push!(x.axes[1], a)
        push!(x.data, @variable(model, lower_bound = 0, upper_bound = 1, base_name = "x[$a]"))
        x.lookup[1].data[a] = length(x.data)

        # Update the objective coefficient of the new column
        set_objective_coefficient(model, x[a], time(data, a))

        # Update constraints
        set_normalized_coefficient(degree[i], x[a], - 1)
        set_normalized_coefficient(degree[j], x[a], 1)
        set_normalized_coefficient(degree_at_most_once[i], x[a], -1)

        selected_blocks::VVi = nodes_blocks[i]
        for block::Vi in filter(block::Vi -> a in δ⁺(A, block), selected_blocks)
            set_normalized_coefficient(block_out_degree[block], x[a], 1)
        end

        for (S::Vi, i_::Int, j_::Int) in subtour_components 
            if a in δ⁺(A′′, S)
                set_normalized_coefficient(subtour[Tuple{Vi, Int, Int}((S, i_, j_))], x[a], 1)
            end

            if i == i_ || i == j_
                set_normalized_coefficient(subtour[Tuple{Vi, Int, Int}((S, i_, j_))], x[a], -1)
            end
        end

    end

    function appendSubtourConstraint(S::Vi, i::Int, j::Int)
        push!(subtour.axes[1], Tuple{Vi, Int, Int}((S, i, j)))
        push!(
              subtour.data, 
              @constraint(
                          model, 
                          sum(a::Arc -> x[a], δ⁺(A′′, S)) - sum(a::Arc -> x[a], δ⁺(A′′, i)) - sum(a::Arc -> x[a], δ⁺(A′′, j)) + 1 >= 0, 
                          base_name = "subtour[$S, $i, $j]"
                         )
             )
        subtour.lookup[1].data[Tuple{Vi, Int, Int}((S, i, j))] = length(subtour.data)
    end

    # generate columns
    curr_iteration::Int = 0
    num_columns_per_it::Int = 100
    max_iterations::Int = 200
    while curr_iteration < max_iterations

        # add subtour cuts
        while true

            optimize!(model)

            #
            subtour_component::Vector{Tuple{Vi, Int, Int}} = getSubtourCuts(V, A′′, model, app)

            #
            isempty(subtour_component) && break

            #
            for (S::Vi, i::Int, j::Int) in subtour_component
                appendSubtourConstraint(S, i, j)

                push!(subtour_components, Tuple{Vi, Int, Int}((S, i, j)))
            end

        end

        # error checking 1
#        if termination_status(model) != MOI.OPTIMAL
        if !in(termination_status(model), [MOI.OPTIMAL, MOI.TIME_LIMIT, MOI.ALMOST_INFEASIBLE]) 
            throw(InvalidStateException("The model could not be solved"))
        end

        # error checking 2
        for a::Arc in A′′
            if abs(JuMP.reduced_cost(x[a]) - getDual(a)) > 1e-5
                println("Error in arc $a ($(value(x[a]))): $(JuMP.reduced_cost(x[a])) != $(getDual(a))")
#                println("$(dual(LowerBoundRef(x[a])))")
#                println("$(dual(UpperBoundRef(x[a])))")
            end
        end

        # filter arcs
        candidate_arcs::Arcs = setdiff(A, A′′)

        # edge case
        isempty(candidate_arcs) && break

        # get best arcs
        arcs_duals::Vector{Tuple{Float64, Arc}} = Vector{Tuple{Float64, Arc}}(map(
                                                                                  a::Arc -> Tuple{Float64, Arc}((getDual(a), a)), 
                                                                                  candidate_arcs
                                                                                 ))

#        filter!((dual, _)::Tuple{Float64, Arc} -> dual <= 0, arcs_duals)
        filter!((dual, _)::Tuple{Float64, Arc} -> dual < 0, arcs_duals)

        sort!(arcs_duals)

        selected_arcs::Arcs = map(
                                  (_, a)::Tuple{Float64, Arc} -> a, 
                                  arcs_duals[begin:min(num_columns_per_it, length(arcs_duals))]
                                 )

        # error checking 3
        for a in selected_arcs
            dual = getDual(a)
#            if dual > 0
            if dual >= 0
                println("==Error")
            end
        end

        # edge case
        if isempty(selected_arcs)
            @debug @sprintf("No available variable with reduced cost")
            break
        end

        # store arc

        # log
        println(@sprintf("Iteration %d", curr_iteration))
        println(@sprintf("\tObjective function: %f", objective_value(model)))
        println(@sprintf("\t# Arcs: %d", length(A′′)))
        println(@sprintf("\t# Selected arcs: %d", length(selected_arcs)))
        
        for selected_arc::Arc in selected_arcs
            if selected_arc in A′′
                println("===Error here")
            end
            push!(A′′, selected_arc)
            appendVariable(selected_arc)
        end

        # increment
        curr_iteration += 1
    end

    # Branch and price

    # creating model with both variables (x and y) relaxed
    set_binary.(x)
    unset_silent(model)

    optimize!(model)

    println("Objective function: ", objective_value(model))
    println("Relative GAP: ", relative_gap(model))

    for a::Arc in A′′
        i, j = a
        if value(x[a]) > 0.5
            println(i, ", ", j)
        end
    end

    @debug "Getting initial relaxation"

#    info["initialLP"] = string(...)

    @debug "Getting initial relaxation with w as integer"

#    info["wLPTime"] = string(@elapsed optimize!(model))
#    info["wLP"] = string(objective_value(model))

    # run
#    info["solverTime"]  = string(@elapsed optimize!(model))
    info["cost"]        = @sprintf("%.2f", objective_value(model))
    info["relativeGAP"] = string(relative_gap(model))
    info["nodeCount"]   = string(node_count(model))

    # retrieve solution
    solution_arcs::Arcs = Arcs(filter(a::Arc -> value(x[a]) > 0.5, A′′))

    tour::Vi = Vi([first(first(solution_arcs))])

    while !isempty(solution_arcs)

        arc_idx::Union{Int, Nothing} = findfirst(a::Arc -> first(a) == last(tour), solution_arcs)

        if arc_idx == nothing 
            error("Desired arc not found")
        end

        a::Arc = solution_arcs[arc_idx]

        i::Int, j::Int = first(a), last(a)

        push!(tour, j)

        deleteat!(solution_arcs, arc_idx)
    end

    solution::SBRPSolution   = SBRPSolution(tour, B)

    # return
    return solution, info

end
