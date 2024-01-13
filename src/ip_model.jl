using Printf
using Combinatorics

include("SparseMaxFlowMinCut.jl")

#=
Add subtour inequalities 
input: 
- model::Model is a Mathematical Programming model
- sets::Set{Pair{Arcs, Arcs}} is a list of list of arcs
=# 
function addSubtourCuts(model::Model, sets::Set{Tuple{Arcs, Arcs}})

    # get vars
    x = model[:x]

    @debug "Adding subtour cuts"

    for (lhs_arcs::Arcs, rhs_arcs::Arcs) in sets
        # add constraint
        @constraint(model, sum(lhs_arc::Arc -> x[lhs_arc], lhs_arcs) >= sum(rhs_arc::Arc -> x[rhs_arc], rhs_arcs) - 1)
    end
end

#=
A Lazy Separation
input: 
- data::SBRPData is a SBRP instance
- Vb::Set{Int} is the set of nodes belonging to some block
- info::Dict{String, String} is the output log relation
- model::Model is a Mathematical Programming model
- cb_data::CPLEX.CallbackContext is the callback context
- context_id::Clong is the context id
=# 
#=
function lazySeparation(
        data::SBRPData, 
        Vb::Si, 
        info::Dict{String, String}, 
        model::Model, 
        cb_data::CPLEX.CallbackContext, 
        context_id::Clong
    )

    @debug "Lazy separation"

    x, y = model[:x], model[:y]

    # preliminary checkings
    context_id != CPX_CALLBACKCONTEXT_CANDIDATE && return

    ispoint_p::Ref{Cint} = Ref{Cint}()

    (CPXcallbackcandidateispoint(cb_data, ispoint_p) != 0 || ispoint_p[] == 0) && return 

    # get values
    CPLEX.load_callback_variable_primal(cb_data, context_id)

    x_val, y_val = callback_value.(Ref(cb_data), x), callback_value.(Ref(cb_data), y)

    # bfs
    A::Arcs = data.D.A
    A′ = filter(a::Arc -> x_val[a] > 0.5, A)
    sets = bfs_sets(A′, Vb, data.depot)

    # get valid components
    #  components = Set{Pair{Si, Int, Vi}}((S, i, block) for S in sets for block in blocks(data.B, S) for i in intersect(block, S) if z_val[i] + y_val[block] - 1 > length(δ⁺(A′, S)) + 0.5) 
    components::Set{Pair{Arcs, Arcs}} = Set{Pair{Arcs, Arcs}}()

    for S in sets 

        lhs = sum(a::Arc -> x_val[a], δ⁺(A′, S))
        Aₛ, A′ₛ = δ⁺(A, S), δ⁺(A′, S) 

        for i::Int in S 

            A′ᵢ = δ⁺(A′, i)

            if i in Vb && sum(a::Arc -> x_val[a], A′ₛ) + 0.5 < sum(a::Arc -> x_val[a], A′ᵢ)

                Aᵢ = δ⁺(A, i)

                push!(components, (Aₛ, Aᵢ))
                MOI.submit(model, MOI.LazyConstraint(cb_data), @build_constraint(sum(a::Arc -> x[a], Aₛ) >= sum(a::Arc -> x[a], Aᵢ))) 

            end
        end
    end

    # add ineqs
    #  [MOI.submit(model, MOI.LazyConstraint(cb_data), @build_constraint(Σ(x[a] for a in δ⁺(A, S)) >= z[i] + y[block] - 1)) for (S, i, block) in components]

    # update info
    info["lazyCuts"] = string(parse(Int, info["lazyCuts"]) + length(components))
end
=#

function addSubtourCut(model::Model, lhs_arcs::Arcs, rhs_arcs::Arcs)

    # get vars
    x = model[:x]

    # add constraint
    @constraint(model, sum(lhs_arc::Arc -> x[lhs_arc], lhs_arcs) >= sum(rhs_arc::Arc -> x[rhs_arc], rhs_arcs) - 1)

end
#=
Get subtour cuts
input: 
- data::SBRPData is the instance
- model::Model is a Mathematical Programming model
- info::Dict{String, String} is the output log relation
output:
- components::Set{Tuple{Arcs, Arcs}} is the set of components to separate.
=# 
function getSubtourCuts(data::SBRPData, model::Model, app::Dict{String, Any}, info::Dict{String, String})::Set{Tuple{Arcs, Arcs}}

    @debug "Obtaining subtour cuts"

    x, z, w = model[:x], model[:z], model[:w]

    # setup
    V::Vi = collect(keys(data.D.V))
    A::Arcs = data.D.A
    B::VVi = data.B

    # outputs
    components::Set{Tuple{Arcs, Arcs}} = Set{Tuple{Arcs, Arcs}}()

    # helpers
    Vₘ = Dict{Int, Int}(map((idx, i)::Tuple{Int, Int} -> i => idx, enumerate(V)))
    Vₘʳ = Dict{Int, Int}(map((idx, i)::Tuple{Int, Int} -> idx => i, enumerate(V)))

    n::Int = length(Vₘ)
    iteration::Int = 1

    # get cuts greedly
    while true
        # store time
        time::Float64 = @elapsed optimize!(model)

        # error checking
        if !in(termination_status(model), [MOI.OPTIMAL, MOI.TIME_LIMIT, MOI.ALMOST_INFEASIBLE]) 
            throw(InvalidStateException("The model could not be solved"))
        end

        # get values
        w_val::Dict{Int, Float64} = Dict{Int, Float64}(map(i::Int -> i => value(w[i]), V))
        z_val::Dict{Int, Float64} = Dict{Int, Float64}(map(i::Int -> i => value(z[i]), V))
        x_val::ArcCostMap = ArcCostMap(map(a::Arc -> a => value(x[a]), A))

        # get subsets
        g = SparseMaxFlowMinCut.ArcFlow[]
        M::Int = 100000
        new_components::Set{Tuple{Arcs, Arcs}} = Set{Tuple{Arcs, Arcs}}()

        # used arcs
        A′::Arcs = filter(a::Arc -> x_val[a] > EPS, A)

        # used nodes
        V′::Vi = filter(i::Int -> z_val[i] > EPS, V)

        # depots
        depots′::Vi = filter(i::Int -> w_val[i] > EPS, V)

        # mounting graph
        for (i::Int, j::Int) in A′
            push!(g, SparseMaxFlowMinCut.ArcFlow(Vₘ[i], Vₘ[j], trunc(floor(x_val[Arc(i, j)], digits=5) * M)))
        end

        # 
        max_violation::Float64 = 0.0

        # get subsets
        for source::Int in depots′

            for target::Int in V′

                # edge case 
                source == target && continue

                # init
                maxFlow::Float64, flows, set = SparseMaxFlowMinCut.find_maxflow_mincut(
                                                                                       SparseMaxFlowMinCut.Graph(n, g), 
                                                                                       Vₘ[source], 
                                                                                       Vₘ[target]
                                                                                      )
                flow::Float64 = maxFlow / M

                # base case: In the same component
                (set[Vₘ[target]] == 1) && continue

                # base case: Condition not met
                flow + 1e-2 >= z_val[source] + z_val[target] - 1 && continue

                # get set
                S::Si = Si(map(i::Int -> Vₘʳ[i], filter(i::Int -> set[i] == 1, 1:n)))

                # base case 3
                #      length(S) <= 1 && continue

                # get components
                Aₛ::Arcs = δ⁺(A, S)
                Aᵢ::Arcs = union(δ⁺(A′, source), δ⁺(A′, target))

                # calculate violation degree
                violation::Float64 = z_val[source] + z_val[target] - 1 - (flow + 1e-2)


                # best improvement case
                if app["subcycle-separation"] == "best" 

                    if max_violation < violation

                        # clean
                        empty!(new_components)

                    else

                        continue

                    end
                end

                # update
                max_violation = max(violation, max_violation)

                # store
                push!(new_components, (Aₛ, Aᵢ))
                addSubtourCut(model, Aₛ, Aᵢ)

                # first improvement case
                app["subcycle-separation"] == "first" && break

            end

            # first improvement case
            (app["subcycle-separation"] == "first" && !isempty(new_components)) && break

        end


        # base case
        isempty(new_components) && break

        # update infos
        #    info["iteration_" * string(iteration) * "_time"], info["iteration_" * string(iteration) * "_cuts"], iteration = time, length(sets′), iteration + 1

        # store components
        union!(components, new_components)

        # add ineqs
        addSubtourCuts(model, components)
    end
    return components
end

#=
Build SBRP model
input: 
- data::SBRPData is a SBRP instance
- app::Dict{String, Any} is the relation of application parameters
=# 
function runCOPCompleteDigraphIPModel(
        data::SBRPData, 
        app::Dict{String, Any},
        maximal_paths::VVi = VVi(),
        time::Function = time,
        blockTime::Function = blockTime
    )::Tuple{SBRPSolution, Dict{String, String}}

    # instance parameters
    depot::Int                 = data.depot
    B::VVi                     = data.B
    A::Arcs                    = data.D.A
    A_set::ArcsSet             = ArcsSet(data.D.A)
    T::Float64                 = data.T
    V::Vi                      = collect(keys(data.D.V))
    profits::Dict{Vi, Float64} = data.profits
    Vb::Si                     = getBlocksNodes(data)
    info::Dict{String, String} = Dict{String, String}("lazyCuts" => "0") 

    # output
    subtour_cuts::Set{Tuple{Arcs, Arcs}} = Set{Pair{Arcs, Arcs}}()

    # new digraph
    A′::Arcs = filter((i, j)::Arc -> j != depot, A)
    V′::Si = setdiff(Si(V), depot)

    P′::Arcs = filter((i, j)::Arc -> i < j && Arc(j, i) in A_set, A)

    nodes_blocks::Dict{Int, VVi} = Dict{Int, VVi}(map(i::Int -> i => filter(block::Vi -> i in block, B), collect(Vb)))

    function createModel(
            relax_x::Bool = false, 
            relax_y::Bool = false,
            relax_z::Bool = false,
            relax_w::Bool = false
        )::Model

        local model

        model::Model = direct_model(CPLEX.Optimizer())
        set_silent(model)
        set_parameters(model, "CPX_PARAM_TILIM" => 3600)
        #=
        set_parameters(model, "CPX_PARAM_PREIND" => 0)
        set_parameters(model, "CPXPARAM_MIP_Strategy_Search" => 1)
        set_parameters(model, "CPXPARAM_Preprocessing_Fill" => 0)
        set_parameters(model, "CPXPARAM_Preprocessing_Dual" => -1)
        set_parameters(model, "CPXPARAM_Preprocessing_Reduce" => 0)
        set_parameters(model, "CPXPARAM_MIP_Cuts_Covers" => -1)
        set_parameters(model, "CPXPARAM_MIP_Cuts_FlowCovers" => -1)
        set_parameters(model, "CPXPARAM_MIP_Cuts_MIRCut" => -1)
        set_parameters(model, "CPXPARAM_MIP_Cuts_ZeroHalfCut" => -1)
        set_parameters(model, "CPXPARAM_MIP_Cuts_Gomory" => -1)
        set_parameters(model, "CPXPARAM_MIP_Cuts_Implied" => -1)
        set_parameters(model, "CPXPARAM_MIP_Cuts_LiftProj" => -1)
        set_parameters(model, "CPXPARAM_MIP_Cuts_Cliques" => -1)
        =#

        #        set_silent(model)
        if relax_x
            @variable(model, x[a in A], lower_bound = 0, upper_bound = 1)
        else
            @variable(model, x[a in A], Bin)
        end
        if relax_y
            @variable(model, y[b in B], lower_bound = 0, upper_bound = 1)
        else
            @variable(model, y[b in B], Bin)
        end
        if relax_z
            @variable(model, z[i in V], lower_bound = 0, upper_bound = 1)
        else
            @variable(model, z[i in V], Bin)
        end
        if relax_w
            @variable(model, w[i in V], lower_bound = 0, upper_bound = 1)
        else
            @variable(model, w[i in V], Bin)
        end

        @objective(model, Max, sum(block::Vi -> data.profits[block] * y[block], B))

        @constraint(model, sum(i::Int -> w[i], V) == 1)
        @constraint(model, lift_w[i::Int in V], w[i] <= z[i])
        @constraint(model, degree[i::Int in V], sum(a::Arc -> x[a], δ⁻(A, i)) == sum(a::Arc -> x[a], δ⁺(A, i)))
        @constraint(model, update_z[i::Int in V], sum(a::Arc -> x[a], δ⁻(A, i)) == z[i])
        @constraint(model, serviced_block[block::Vi in B], sum(a::Arc -> x[a], δ⁺(A, block)) >= y[block])
        @constraint(model, sum(a::Arc -> time(data, a) * x[a], A) <= T - sum(block::Vi -> y[block] * blockTime(data, block), B))

        # improvements
        @constraint(model, block3[block::Vi in B], y[block] - sum(x[a] for i in block for a in δ⁺(A, i) if length(nodes_blocks[i]) == 1) >= 0)
        @constraint(model, subcycle_size_two[a::Arc in P′], x[a] + x[reverse(a)] <= 1)
        
        # fixed depot
        if !app["unfixed-depot"]
            @constraint(model, z[depot] == 1)
            @constraint(model, w[depot] == 1)
        end

        # MTZs
        if app["arcs-mtz"] # Arcs MTZ
            if app["unfixed-depot"]
                error("Invalid flag 'unfixed-depot' when using 'arcs-mtz'")
            else
                @variable(model, t[a::Arc in A], lower_bound = 0, upper_bound = T)
                @constraint(model, sum(a::Arc -> t[a], δ⁺(A, depot)) == 0.0)
                @constraint(model, mtz[i::Int in V′], sum(a::Arc -> t[a], δ⁺(A, i)) == sum(a::Arc -> t[a], δ⁻(A, i)) + sum(a::Arc -> x[a] * time(data, a), δ⁺(A, i)))
                @constraint(model, ub[i::Int in V′], sum(a::Arc -> t[a], δ⁻(A, i)) <= T - sum(block::Vi -> y[block] * blockTime(data, block), B))
                @constraint(model, ub1[a::Arc in A], t[a] <= x[a] * T)
                @constraint(model, ub2[i::Int in V′], sum(a::Arc -> t[a], δ⁻(A, i)) <= z[i] * T)
            end
        else # Nodes MTZ
            @variable(model, t[i::Int in V], lower_bound = 0, upper_bound = T)
            if app["unfixed-depot"]
                @constraint(model, mtz[a::Arc in filter(a::Arc -> reverse(a) in A_set, A)], t[last(a)] >= t[first(a)] + x[a] * time(data, a) - (1 - x[a]) * T - x[reverse(a)] * time(data, reverse(a)) - w[last(a)] * T)
                @constraint(model, mtz1[a::Arc in filter(a::Arc -> !in(reverse(a), A_set), A)], t[last(a)] >= t[first(a)] + x[a] * time(data, a) - (1 - x[a]) * T - w[last(a)] * T)
                @constraint(model, ub1[i::Int in V], t[i] <= T - sum(block::Vi -> y[block] * blockTime(data, block), B))
                @constraint(model, ub2[i::Int in V], t[i] <= z[i] * T)
                @constraint(model, ub3[i::Int in V], t[i] <= (1 - w[i]) * T)
            else
                @constraint(model, t[depot] == 0.0)
                @constraint(model, mtz[a::Arc in A′], t[last(a)] >= t[first(a)] + x[a] * time(data, a) - (1 - x[a]) * T - x[reverse(a)] * time(data, reverse(a)))
                @constraint(model, ub1[i::Int in V′], t[i] <= T - sum(block::Vi -> y[block] * blockTime(data, block), B))
                @constraint(model, ub2[i::Int in V′], t[i] <= z[i] * T)
            end
        end
        # subtour cuts
#        addSubtourCuts(model, subtour_cuts) 
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
        model[:x], model[:y], model[:z], model[:w], model[:t] = x, y, z, w, t

        return model
    end

    # creating model with both variables (x and y) relaxed
    model::Model = createModel(true, true)

    # getting intersection cuts
    if app["intersection-cuts"]
        @debug "Getting intersection cuts"

        # edge case
        if isempty(maximal_paths)
            error("The maximal paths list is empty")
        end

        # add inequalities
        z = model[:z]
        for maximal_path::Vi in maximal_paths
            @constraint(model, sum(i::Int -> z[i], maximal_path) <= 1)
        end
    end

    # getting initial relaxation with both variables <x and y> relaxed
    @debug "Getting initial relaxation"
    optimize!(model)
    info["initialLP"] = string(objective_value(model))

    # getting initial relaxation with only x relaxed (y integer)
    @debug "Getting initial relaxation with y as integer"
    model = createModel(true)

    info["yLPTime"] = string(@elapsed optimize!(model))
    info["yLP"] = string(objective_value(model))

    # get max-flow cuts with x and y relaxed or integer
    if app["subcycle-separation"] != "none"

        model = createModel(true, !app["y-integer"], !app["z-integer"], !app["w-integer"])

        optimize!(model)

        info["maxFlowCutsTime"] = string(@elapsed subtour_cuts = getSubtourCuts(data, model, app, info))
        info["maxFlowCuts"] = string(length(subtour_cuts))
        info["maxFlowLP"] = string(objective_value(model))
    end

    # integer model
    model = createModel()
    #  MOI.set(model, MOI.NumberOfThreads(), 1) # thread numbers
    #  MOI.set(model, CPLEX.CallbackFunction(), (cb_data, context_id) -> lazy_separation(data, Vb, info, model, cb_data, context_id)) # lazy callback

    # run
    info["solverTime"]  = string(@elapsed optimize!(model))
    info["cost"]        = @sprintf("%.2f", objective_value(model))
    info["relativeGAP"] = string(relative_gap(model))
    info["nodeCount"]   = string(node_count(model))

    # retrieve solution
    x, y, z, w, t = model[:x], model[:y], model[:z], model[:w], model[:t]

    solution_arcs::Arcs  = Arcs(filter(a::Arc -> value(x[a]) > 0.5, A))
    solution_nodes::Vi   = Vi(filter(i::Int -> value(z[i]) > 0.5, V))
    solution_blocks::VVi = VVi(filter(block::Vi -> value(y[block]) > 0.5, B))
    chosen_depot::Int    = depot

    if app["unfixed-depot"]
        chosen_depot = first(filter(i::Int -> value(w[i]) > 0.5, V))
    end

    tour::Vi = Vi([chosen_depot])

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
end

#=
Build CSP model
input: 
- data::SBRPData is a CSP instance
- app::Dict{String, Any} is the relation of application parameters
=# 
function runCSPCompleteDigraphIPModel(
        data::SBRPData, 
        app::Dict{String, Any}, 
        maximal_paths::VVi = VVi(),
        time::Function = time,
        blockTime::Function = blockTime
    )::Tuple{SBRPSolution, Dict{String, String}}

    # instance parameters
    depot::Int                 = data.depot
    B::VVi                     = data.B
    A::Arcs                    = data.D.A
    A_set::ArcsSet             = ArcsSet(data.D.A)
    T::Float64                 = data.T
    V::Vi                      = collect(keys(data.D.V))
    profits::Dict{Vi, Float64}   = data.profits
    distance::Dict{Arc, Float64} = data.D.distance
    info::Dict{String, String} = Dict{String, String}("lazyCuts" => "0") 

    # output
    subtour_cuts::Set{Tuple{Arcs, Arcs}} = Set{Tuple{Arcs, Arcs}}()

    # new digraph
    A′::Arcs = filter((i, j)::Arc -> j != depot, A)
    V′::Si = setdiff(Si(V), depot)

    P′::Arcs = filter((i, j)::Arc -> i < j && Arc(j, i) in A_set, A)

    nodes_blocks::Dict{Int, VVi} = Dict{Int, VVi}(map(i::Int -> i => filter(block::Vi -> i in block, B), V))

    function createModel(
            relax_x::Bool = false, 
            relax_z::Bool = false,
            relax_w::Bool = false
        )::Model

        local model

        model::Model = direct_model(CPLEX.Optimizer())
#        set_silent(model)
        set_parameters(model, "CPX_PARAM_TILIM" => 3600)
        #=
        set_parameters(model, "CPX_PARAM_PREIND" => 0)
        set_parameters(model, "CPXPARAM_MIP_Strategy_Search" => 1)
        set_parameters(model, "CPXPARAM_Preprocessing_Fill" => 0)
        set_parameters(model, "CPXPARAM_Preprocessing_Dual" => -1)
        set_parameters(model, "CPXPARAM_Preprocessing_Reduce" => 0)
        set_parameters(model, "CPXPARAM_MIP_Cuts_FlowCovers" => -1)
        set_parameters(model, "CPXPARAM_MIP_Cuts_MIRCut" => -1)
        set_parameters(model, "CPXPARAM_MIP_Cuts_ZeroHalfCut" => -1)
        set_parameters(model, "CPXPARAM_MIP_Cuts_Gomory" => -1)
        set_parameters(model, "CPXPARAM_MIP_Cuts_Implied" => -1)
        set_parameters(model, "CPXPARAM_MIP_Cuts_LiftProj" => -1)
        set_parameters(model, "CPXPARAM_MIP_Cuts_Cliques" => -1)
        =#

        if relax_x
            @variable(model, x[a in A], lower_bound = 0, upper_bound = 1)
        else
            @variable(model, x[a in A], Bin)
        end
        if relax_z
            @variable(model, z[i in V], lower_bound = 0, upper_bound = 1)
        else
            @variable(model, z[i in V], Bin)
        end
        if relax_w
            @variable(model, w[i in V], lower_bound = 0, upper_bound = 1)
        else
            @variable(model, w[i in V], Bin)
        end

        @objective(model, Min, sum(a::Arc -> time(data, a) * x[a], A))

        @constraint(model, sum(i::Int -> w[i], V) == 1)
        @constraint(model, lift_w[i::Int in V], w[i] <= z[i])
        @constraint(model, degree[i::Int in V], sum(a::Arc -> x[a], δ⁻(A, i)) == sum(a::Arc -> x[a], δ⁺(A, i)))
        @constraint(model, degree_at_most_once[i::Int in V], sum(a::Arc -> x[a], δ⁺(A, i)) == z[i])
        @constraint(model, block_out_degree[block::Vi in B], sum(a::Arc -> x[a], δ⁺(A, block)) >= 1)
#        @constraint(model, block_out_degree[block::Vi in B], sum(i::Int -> z[i], block) >= 1)
        @constraint(model, no_two_subcycle[a::Arc in P′], x[a] + x[reverse(a)] <= 1)

        # fixed depot
        if !app["unfixed-depot"]
            @constraint(model, z[depot] == 1)
            @constraint(model, w[depot] == 1)
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
                @constraint(model, ub1[i::Int in V′],  sum(a::Arc -> t[a], δ⁻(A, i)) <= z[i] * T)
            end
        else # Nodes MTZ
            @variable(model, t[i::Int in V], lower_bound = 0, upper_bound = T)
            if app["unfixed-depot"]
                @constraint(model, mtz[a::Arc in filter(a::Arc -> reverse(a) in A_set, A)], t[last(a)] >= t[first(a)] + x[a] * time(data, a) - (1 - x[a]) * T - x[reverse(a)] * time(data, reverse(a)) - w[last(a)] * T)
                @constraint(model, mtz1[a::Arc in filter(a::Arc -> !in(reverse(a), A_set), A)], t[last(a)] >= t[first(a)] + x[a] * time(data, a) - (1 - x[a]) * T - w[last(a)] * T)
                @constraint(model, ub1[i::Int in V], t[i] <= z[i] * T)
                @constraint(model, ub2[i::Int in V], t[i] <= (1 - w[i]) * T)
            else
                @constraint(model, t[depot] == 0.0)
                @constraint(model, mtz[a::Arc in filter(a::Arc -> reverse(a) in A_set, A′)], t[last(a)] >= t[first(a)] + x[a] * time(data, a) - (1 - x[a]) * T - x[reverse(a)] * time(data, reverse(a)) - w[last(a)] * T)
                @constraint(model, mtz1[a::Arc in filter(a::Arc -> !in(reverse(a), A_set), A′)], t[last(a)] >= t[first(a)] + x[a] * time(data, a) - (1 - x[a]) * T - w[last(a)] * T)
                @constraint(model, ub1[i::Int in V′], t[i] <= z[i] * T)
            end
        end

        # subtour cuts
        addSubtourCuts(model, subtour_cuts) 

        # lb
        if app["lb"] != nothing 
            lb::Int = parse(Int, app["lb"])
            @constraint(model, sum(a::Arc -> distance[a] * x[a], A) >= lb)
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
        model[:x], model[:z], model[:w] = x, z, w

        return model
    end

    # creating model with both variables (x and y) relaxed
    model::Model = createModel(true, true)

    # getting intersection cuts
    if app["intersection-cuts"]
        @debug "Getting intersection cuts"

        # edge case
        if isempty(maximal_paths)
            error("The maximal paths list is empty")
        end
        
        # add inequalities
        z = model[:z]
        for maximal_path::Vi in maximal_paths
            @constraint(model, sum(i::Int -> z[i], maximal_path) <= 1)
        end
    end

    # getting initial relaxation with both variables <x and y> relaxed
    @debug "Getting initial relaxation"
    optimize!(model)
    info["initialLP"] = string(objective_value(model))

    # getting initial relaxation with only x relaxed (y integer)
    @debug "Getting initial relaxation with y as integer"
    model = createModel(true)

    info["zLPTime"] = string(@elapsed optimize!(model))
    info["zLP"] = string(objective_value(model))

    # get max-flow cuts with x and y relaxed or integer
    if app["subcycle-separation"] != "none"

        model = createModel(true, !app["z-integer"], !app["w-integer"])
        
        optimize!(model)

        info["maxFlowCutsTime"] = string(@elapsed subtour_cuts = getSubtourCuts(data, model, app, info))
        info["maxFlowCuts"] = string(length(subtour_cuts))
        info["maxFlowLP"] = string(objective_value(model))
    end

    # integer model
    model = createModel()
    #  MOI.set(model, MOI.NumberOfThreads(), 1) # thread numbers
    #  MOI.set(model, CPLEX.CallbackFunction(), (cb_data, context_id) -> lazy_separation(data, Vb, info, model, cb_data, context_id)) # lazy callback

    # run
    info["solverTime"]  = string(@elapsed optimize!(model))
    info["cost"]        = @sprintf("%.2f", objective_value(model))
    info["relativeGAP"] = string(relative_gap(model))
    info["nodeCount"]   = string(node_count(model))

    # retrieve solution
    x, z, w = model[:x], model[:z], model[:w]

    solution_arcs::Arcs = Arcs(filter(a::Arc -> value(x[a]) > 0.5, A))
    visited_nodes::Vi   = filter(i::Int -> value(z[i]) > 0.5, V)
    chosen_depot::Int = first(filter(i::Int -> value(w[i]) > 0.5, V))

    tour::Vi = Vi([chosen_depot])

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
