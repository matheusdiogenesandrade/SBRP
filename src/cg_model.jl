using Printf

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
function runCSPColumnGenerationModel(data::SBRPData, app::Dict{String, Any})::Tuple{SBRPSolution, Dict{String, String}}

    # instance parameters
    depot::Int                 = data.depot
    B::VVi                     = data.B
    A::Arcs                    = data.D.A
    A_set::ArcsSet             = ArcsSet(data.D.A)
    T::Float64                 = data.T
    V::Dict{Int, Vertex}       = data.D.V
    Vb::Si                     = getBlocksNodes(data)
    nodes_blocks::Dict{Int, VVi} = Dict{Int, VVi}(map(i::Int -> i => filter(block::Vi -> i in block, B), collect(Vb)))

    # nodes indexes
    V_seq::Vi = collect(keys(V))
    V_idx::Dict{Int, Int} = Dict{Int, Int}(i => idx for (idx::Int, i::Int) in enumerate(V_seq))

    # util
    P′::Arcs = filter((i, j)::Arc -> i < j && Arc(j, i) in A_set, A)

    # output
    info::Dict{String, String} = Dict{String, String}("lazyCuts" => "0") 

    # MTZ data
#    A′::Arcs = filter((i, j)::Arc -> j != depot, A)
#    V′::Si = setdiff(keys(V), depot)
#    P′::Arcs = filter((i, j)::Arc -> i < j, χ(V′))

    # Initial arcs
    m::Int = length(B)

    A′′::Arcs = []

    if !app["unfixed-depot" ]

        for block::Vi in B
            if depot != first(block) 
                push!(A′′, Arc(depot, first(block)))
                push!(A′′, Arc(first(block), depot))
            end
        end

    end

    for (block1::Vi, block2::Vi) in zip(B[1:m - 1], B[2:m])
        first(block1) != first(block2) && push!(A′′, Arc(first(block1), first(block2)))
    end
    first(last(B)) != first(first(B)) && push!(A′′, Arc(first(last(B)), first(first(B))))

    V′′::Vi = collect(filter(i::Int -> any(a::Arc -> a[1] == i || a[2] == i, A′′), keys(V)))

    # helpers

    function getDual(a::Arc, model::Model)::Float64
        # registries
        # variables
        model[:x], model[:y], model[:t] = x, y, t

        # constraints
        model[:depot_choose] = depot_choose
        model[:lift_w] = lift_w
        model[:degree] = degree
        model[:degree_at_most_once] = degree_at_most_once
        model[:block_out_degree] = block_out_degree
        if !app["unfixed-depot"]
            model[:fixed_depot] = fixed_depot
        end
        if app["ub"] != nothing 
            model[:upper_bound_obj_func] = upper_bound_obj_func
        end
        if app["lb"] != nothing 
            model[:lower_bound_obj_func] = lower_bound_obj_func
        end
        # mtz
        model[:mtz] = mtz
        model[:mtz1] = mtz1
        model[:mtz_ub1] = mtz_ub1
        model[:mtz_ub2] = mtz_ub2

        ##############################
        # params
        i::Int, j::Int = a
        selected_blocks::VVi = i == depot ? [] : nodes_blocks[i]

        i_key::JuMP.Containers.DenseAxisArrayKey{Tuple{Int64}} = JuMP.Containers.DenseAxisArrayKey{Tuple{Int64}}((i,))
        j_key::JuMP.Containers.DenseAxisArrayKey{Tuple{Int64}} = JuMP.Containers.DenseAxisArrayKey{Tuple{Int64}}((j,))

        degree_i_dual::Float64              = i_key in keys(model[:degree]) ? - dual(model[:degree][i_key]) : 0.0
        degree_j_dual::Float64              = j_key in keys(model[:degree]) ?   dual(model[:degree][j_key]) : 0.0
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

    function createModel()::Model

        local model

        model::Model = direct_model(CPLEX.Optimizer())
#        set_silent(model)
#        set_parameters(model, "CPX_PARAM_TILIM" => 3600)
        #=
        set_parameters(model, "CPX_PARAM_PREIND" => 0)
        set_parameters(model, "CPXPARAM_MIP_Strategy_Search" => 1)
        set_parameters(model, "CPXPARAM_Preprocessing_Fill" => 0)
        set_parameters(model, "CPXPARAM_Preprocessing_Dual" => -1)
        set_parameters(model, "CPXPARAM_Preprocessing_Reduce" => 0)
        =#
        set_parameters(model, "CPXPARAM_MIP_Cuts_FlowCovers" => -1)
        set_parameters(model, "CPXPARAM_MIP_Cuts_MIRCut" => -1)
        set_parameters(model, "CPXPARAM_MIP_Cuts_ZeroHalfCut" => -1)
        set_parameters(model, "CPXPARAM_MIP_Cuts_Gomory" => -1)
        set_parameters(model, "CPXPARAM_MIP_Cuts_Implied" => -1)
        set_parameters(model, "CPXPARAM_MIP_Cuts_LiftProj" => -1)
        set_parameters(model, "CPXPARAM_MIP_Cuts_Cliques" => -1)

        @variable(model, x[a in A], lower_bound = 0, upper_bound = 1)
#        @variable(model, z[i in V], lower_bound = 0, upper_bound = 1)
        @variable(model, w[i in V], lower_bound = 0, upper_bound = 1)

        @objective(model, Min, sum(a::Arc -> time(data, a) * x[a], A))

        @constraint(model, depot_choose, sum(i::Int -> w[i], V) == 1)
#        @constraint(model, lift_w[i::Int in V], w[i] <= z[i])
        @constraint(model, lift_w[i::Int in V], w[i] <= sum(a::Arc -> x[a], δ⁺(A, i)))
        @constraint(model, degree[i::Int in V], sum(a::Arc -> x[a], δ⁻(A, i)) == sum(a::Arc -> x[a], δ⁺(A, i)))
        @constraint(model, degree_at_most_once[i::Int in V], sum(a::Arc -> x[a], δ⁺(A, i)) <= 1)
#        @constraint(model, degree_at_most_once[i::Int in V], sum(a::Arc -> x[a], δ⁺(A, i)) == z[i])
        @constraint(model, block_out_degree[block::Vi in B], sum(a::Arc -> x[a], δ⁺(A, block)) >= 1)
#        @constraint(model, block_out_degree[block::Vi in B], sum(i::Int -> z[i], block) >= 1)
#        @constraint(model, no_two_subcycle[a::Arc in P′], x[a] + x[reverse(a)] <= 1)

        # fixed depot
        if !app["unfixed-depot"]
#            @constraint(model, z[depot] == 1)
            @constraint(model, fixed_depot, w[depot] == 1)
        end

        # MTZs
        if app["arcs-mtz"] # Arcs MTZ
            if app["unfixed-depot"]
                error("Invalid flag 'unfixed-depot' when using 'arcs-mtz'")
            else
                @variable(model, t[a::Arc in A], lower_bound = 0, upper_bound = T)
                @constraint(model, ub[a::Arc in A], t[a] <= x[a] * T)
                @constraint(model, sum(a::Arc -> t[a], δ⁺(A, depot)) == 0.0)
                @constraint(model, mtz[i::Int in V′],  sum(a::Arc -> t[a], δ⁺(A, i)) == sum(a::Arc -> t[a], δ⁻(A, i)) + sum(a::Arc -> x[a] * time(data, a), δ⁺(A, i)))
                @constraint(model, ub1[i::Int in V′],  sum(a::Arc -> t[a], δ⁻(A, i)) <= sum(a::Arc -> x[a], δ⁻(A, i)) * T)
            end
        else # Nodes MTZ
            @variable(model, t[i::Int in V], lower_bound = 0, upper_bound = T)
            if app["unfixed-depot"]
                @constraint(model, mtz[a::Arc in filter(a::Arc -> reverse(a) in A_set, A)], t[last(a)] >= t[first(a)] + x[a] * time(data, a) - (1 - x[a]) * T - x[reverse(a)] * time(data, reverse(a)) - w[last(a)] * T)
                @constraint(model, mtz1[a::Arc in filter(a::Arc -> !in(reverse(a), A_set), A)], t[last(a)] >= t[first(a)] + x[a] * time(data, a) - (1 - x[a]) * T - w[last(a)] * T)
                @constraint(model, mtz_ub1[i::Int in V], t[i] <= sum(a::Arc -> x[a], δ⁻(A, i)) * T)
                @constraint(model, mtz_ub2[i::Int in V], t[i] <= (1 - w[i]) * T)
            else
                @constraint(model, mtz[a::Arc in filter(a::Arc -> reverse(a) in A_set, A′)], t[last(a)] >= t[first(a)] + x[a] * time(data, a) - (1 - x[a]) * T - x[reverse(a)] * time(data, reverse(a)) - w[last(a)] * T)
                @constraint(model, mtz1[a::Arc in filter(a::Arc -> !in(reverse(a), A_set), A′)], t[last(a)] >= t[first(a)] + x[a] * time(data, a) - (1 - x[a]) * T - w[last(a)] * T)
                @constraint(model, mtz_ub1[i::Int in V′], t[i] <= sum(a::Arc -> x[a], δ⁻(A, i)) * T)
                @constraint(model, mtz_ub2[i::Int in V], t[i] <= (1 - w[i]) * T)
            end
        end

        # subtour cuts
#        addSubtourCuts(model, subtour_cuts) 

        # lb
        if app["lb"] != nothing 
            lb::Int = parse(Int, app["lb"])
            @constraint(model, lower_bound_obj_func, sum(a::Arc -> distance[a] * x[a], A) >= lb)
        end

        # ub
        if app["ub"] != nothing 
            ub::Int = parse(Int, app["ub"])
            @constraint(model, upper_bound_obj_func, sum(a::Arc -> distance[a] * x[a], A) <= ub)
        end

        # registries
        # variables
        model[:x], model[:y], model[:t] = x, y, t

        # constraints
        model[:depot_choose] = depot_choose
        model[:lift_w] = lift_w
        model[:degree] = degree
        model[:degree_at_most_once] = degree_at_most_once
        model[:block_out_degree] = block_out_degree
        if !app["unfixed-depot"]
            model[:fixed_depot] = fixed_depot
        end
        if app["ub"] != nothing 
            model[:upper_bound_obj_func] = upper_bound_obj_func
        end
        if app["lb"] != nothing 
            model[:lower_bound_obj_func] = lower_bound_obj_func
        end
        # mtz
        model[:mtz] = mtz
        model[:mtz1] = mtz1
        model[:mtz_ub1] = mtz_ub1
        model[:mtz_ub2] = mtz_ub2

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
