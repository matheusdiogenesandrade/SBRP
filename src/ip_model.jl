using Printf

include("SparseMaxFlowMinCut.jl")

intersection_cuts1::Vector{Arcs} = Vector{Arcs}()
intersection_cuts2::Vector{Arcs} = Vector{Arcs}()
subtour_cuts::Set{Pair{Arcs, Arcs}} = Set{Pair{Arcs, Arcs}}()

#=
Add an equality cut
input: 
- model::Model is a Mathematical Programming model
- arcs::Vector{Arcs} is a list of list of arcs
=# 
function addIntersectionCuts1(model::Model, arcs::Vector{Arcs})

    @debug "Adding intersection cuts 1"
    addCuts(model, 
            Vector{Pair{AffExpr, AffExpr}}(map(
                                               A::Arcs -> Pair{AffExpr, AffExpr}(sum(a::Arc -> model[:x][a], A), AffExpr(0)), 
                                               arcs
                                              )))

    addCuts(model, 
            Vector{Pair{AffExpr, AffExpr}}(map(
                                               A::Arcs -> Pair{AffExpr, AffExpr}(AffExpr(0), sum(a::Arc -> model[:x][a], A)), 
                                               arcs
                                              )))
end

#=
Add an inequality cut
input: 
- model::Model is a Mathematical Programming model
- arcs::Vector{Arcs} is a list of list of arcs
=# 
function addIntersectionCuts2(model::Model, arcs::Vector{Arcs})

    @debug "Adding intersection cuts 2"

    addCuts(model, 
            Vector{Pair{AffExpr, AffExpr}}( map(
#                                               A::Arcs -> Pair{AffExpr, AffExpr}(AffExpr(2), sum(a::Arc -> model[:x][a], A)), 
                                               A::Arcs -> Pair{AffExpr, AffExpr}(AffExpr(1), sum(a::Arc -> model[:x][a], A)), 
                                               arcs
                                              )))
end

#=
Add subtour inequalities 
input: 
- model::Model is a Mathematical Programming model
- sets::Set{Pair{Arcs, Arcs}} is a list of list of arcs
=# 
function addSubtourCuts(model::Model, sets::Set{Pair{Arcs, Arcs}})

    @debug "Adding subtour cuts"

    addCuts(model, 
            Vector{Pair{AffExpr, AffExpr}}(map(
                (A1, A2)::Pair{Arcs, Arcs} -> Pair{AffExpr, AffExpr}(sum(a::Arc -> model[:x][a], A1), sum(a::Arc -> model[:x][a] for a in A2)), 
                collect(sets)
               )))
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

#=
Get subtour cuts
input: 
- data::SBRPData is a SBRP instance
- model::Model is a Mathematical Programming model
- info::Dict{String, String} is the output log relation
output:
- components::Set{Pair{Arcs, Arcs}} is the set of components to separate.
=# 
function getSubtourCuts(data::SBRPData, model::Model, info::Dict{String, String})::Set{Pair{Arcs, Arcs}}

    @debug "Obtaining subtour cuts"

    x, y = model[:x], model[:y]

    # setup
    A::Arcs = data.D.A
    B::VVi = data.B
    depot::Int = data.depot
    components::Set{Pair{Arcs, Arcs}} = Set{Pair{Arcs, Arcs}}()
    Vb::Si = getBlocksNodes(data) 

    Vₘ = Dict{Int, Int}(map((idx, i)::Tuple{Int, Int} -> i => idx, enumerate(keys(data.D.V))))
    Vₘʳ = Dict{Int, Int}(map((idx, i)::Tuple{Int, Int} -> idx => i, enumerate(keys(data.D.V))))

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
        y_val::Dict{Vi, Float64} = Dict{Vi, Float64}(map(block::Vi -> block => value(y[block]), B))
        x_val::ArcCostMap = ArcCostMap(map(a::Arc -> a => value(x[a]), A))

        # get subsets
        g = SparseMaxFlowMinCut.ArcFlow[]
        M::Int = 100000
        newComponents::Set{Pair{Arcs, Arcs}} = Set{Pair{Arcs, Arcs}}()

        A′::Arcs = filter(a::Arc -> x_val[a] > EPS, A)

        # mounting graph
        for (i::Int, j::Int) in A′
            push!(g, SparseMaxFlowMinCut.ArcFlow(Vₘ[i], Vₘ[j], trunc(floor(x_val[Arc(i, j)], digits=5) * M)))
        end

        # get subsets
        for source::Int in Vb 

            # init
            maxFlow::Float64, flows, set = SparseMaxFlowMinCut.find_maxflow_mincut(SparseMaxFlowMinCut.Graph(n, g), Vₘ[source], Vₘ[depot])
            flow::Float64 = maxFlow / M

            # base case 1
            (flow >= 1 - 1e-2 || set[Vₘ[depot]] == 1) && continue

            # base case 2
            A′ᵢ = δ⁺(A′, source)

            rhs::Float64 = sum(a::Arc -> x_val[a], A′ᵢ, init = 0)

            # get set
            S::Si = Si(map(i::Int -> Vₘʳ[i], filter(i::Int -> set[i] == 1, 1:n)))

            flow + 1e2 >= rhs && continue

            # base case 3
            #      length(S) <= 1 && continue

            # get components
            #      [push!(newComponents, (S, i, block)) for block in blocks(B, S) for i in intersect(block, S) if z_val[i] + y_val[block] - 1 > flow + 1e-2]
            Aₛ::Si = δ⁺(A, S)
            push!(newComponents, (Aₛ, Aᵢ))
            addCuts(model, Vector{Pair{AffExpr, AffExpr}}([Pair{AffExpr, AffExpr}(sum(a::Arc -> x[a], Aₛ), sum(a::Arc -> x[a], Aᵢ))]))

        end

        # base case
        isempty(newComponents) && break
        #    flush_println(length(sets′))
        #    [flush_println(S, i, block) for (S, i, block) in sets′]

        # update infos
        #    info["iteration_" * string(iteration) * "_time"], info["iteration_" * string(iteration) * "_cuts"], iteration = time, length(sets′), iteration + 1

        # store components
        union!(components, newComponents)

        # add ineqs
        addSubtourCuts(model, components)
    end
    return components
end

#=
Get intersection cuts
input: 
- data::SBRPData is a SBRP instance
=# 
function getIntersectionCuts(data::SBRPData)::Tuple{Vector{Arcs}, Vector{Arcs}}

    @debug "Obtaining intersection cuts"

    # setup
    B::VVi = data.B 
    A::Arcs = data.D.A
    cuts1::Vector{Arcs} = Vector{Arcs}()
    cuts2::Vector{Arcs} = Vector{Arcs}()
    cliques::Vector{VVi} = Vector{VVi}()

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
        B′::VVi = filter(block -> value(z[block]) > 0.5, B)

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

    # get nodes of each block
    Vb::Si = getBlocksNodes(data)
    nodes_blocks::Dict{Int, VVi} = Dict{Int, VVi}(map(i::Int -> i => filter(block::Vi -> i in block, B), Vi(collect(Vb))))

    # get cuts
    for clique::VVi in cliques
        # get intersection
        intersection::Vi = ∩(clique...)

        ######## intersection cuts 1 ########

#        isolated_intersection::Vi = setdiff(intersection, ∪(i for block in setdiff(B, clique) for i in block))
        isolated_intersection::Vi = setdiff(intersection, reduce(vcat, setdiff(B, clique)))
        
        #    isolated_intersection = sediff(intersection, ∪(block... for block ∈ setdiff(B, clique)))
        length(isolated_intersection) > 1 && push!(cuts1, Arcs(filter((i, j)::Arc -> i != j, χ(isolated_intersection))))

        @debug "Isolated clique $isolated_intersection"

        ######## intersection cuts 2 ########

        @debug "Clique length $(length(clique))"

        for block::Vi in clique
            # get arcs incident to nodes that belong exclusively to the block 
            #      independent_arcs = [a for (i, blocks) in nodes_blocks for a in δ⁺(A, i) if ∧(length(blocks) == 1, block ∈ blocks)]
            covered_arcs::Arcs = reduce(
                                        vcat, 
                                        map(
                                            i::Int -> δ⁺(A, i), 
                                            filter(i::Int -> length(nodes_blocks[i]) == 1, block)
                                           ), 
                                        init = Arcs()
                                       )

            # edge case
            #      isempty(independent_arcs) && continue
            isempty(covered_arcs) && continue

            # store
            #      push!(cuts2, Arcs(∪([a for i ∈ intersection for a ∈ δ⁺(A, i)], independent_arcs)))
            push!(cuts2, Arcs(covered_arcs))
        end

    end

    return cuts1, cuts2
end

#=
Build SBRP model
input: 
- data::SBRPData is a SBRP instance
- app::Dict{String, Any} is the relation of application parameters
=# 
function runCompleteDigraphIPModel(data::SBRPData, app::Dict{String, Any})::Tuple{SBRPSolution, Dict{String, String}}

    # instance parameters
    depot::Int                 = data.depot
    B::VVi                     = data.B
    A::Arcs                    = data.D.A
    T::Float64                 = data.T
    V::Dict{Int, Vertex}       = data.D.V
    profits::Dict{Vi, Float64} = data.profits
    Vb::Si                     = getBlocksNodes(data)
    info::Dict{String, String} = Dict{String, String}("lazyCuts" => "0") 

    # new digraph
    A′::Arcs = filter((i, j)::Arc -> j != depot, A)
    V′::Si = setdiff(keys(V), depot)
    
    P′::Arcs = filter((i, j)::Arc -> i < j, χ(V′))

    nodes_blocks::Dict{Int, VVi} = Dict{Int, VVi}(map(i::Int -> i => filter(block::Vi -> i in block, B), collect(Vb)))

    function createModel(relax_x::Bool = false, relax_y::Bool = false)::Model

        global intersection_cuts1
        global intersection_cuts2
        global subtour_cuts
        
        local model

        model::Model = direct_model(CPLEX.Optimizer())
        set_parameters(model, "CPX_PARAM_TILIM" => 3600)

        #    set_silent(model)
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


        @objective(model, Max, sum(block::Vi -> data.profits[block] * y[block], B))

        @constraint(model, degree[i::Int in keys(V)], sum(a::Arc -> x[a], δ⁻(A, i)) == sum(a::Arc -> x[a], δ⁺(A, i)))
        @constraint(model, degree_at_most_once[i::Int in keys(V)], sum(a::Arc -> x[a], δ⁺(A, i)) <= 1)
        @constraint(model, sum(a::Arc -> x[a], δ⁺(A, depot)) == 1)
        @constraint(model, block1[block::Vi in B], sum(a::Arc -> x[a], δ⁺(A, block)) >= y[block])
        @constraint(model, sum(a::Arc -> time(data, a) * x[a], A) <= T - sum(block::Vi -> y[block] * blockTime(data, block), B))

        # improvements
        @constraint(model, block3[block::Vi in B], y[block] - sum(x[a] in δ⁺(A, i) for i in block if length(nodes_blocks[i]) == 1) >= 0)
#        @constraint(model, block4[block::Vi in B], sum((i, j)::Arc -> x[(i, j)], filter((i, j)::Arc -> i in block && j in block && nodes_blocks[i] == nodes_blocks[j] && length(nodes_blocks[j]) == 1, A)) == 0)
        @constraint(model, block4[block::Vi in B], sum(a::Arc -> x[a], filter((i, j)::Arc -> i in block && j in block && (length(nodes_blocks[i]) == 1 || length(nodes_blocks[j]) == 1), A), init = 0.0) == 0)
        @constraint(model, subcycle_size_two[a::Arc in P′], x[a] + x[reverse(a)] <= 1)

        # MTZs
        if app["arcs-mtz"] # Arcs MTZ
            @variable(model, t[a::Arc in A], lower_bound = 0, upper_bound = T)
            @constraint(model, sum(a::Arc -> t[a], δ⁺(A, depot)) == 0.0)
            @constraint(model, mtz[i::Int in V′], sum(a::Arc -> t[a], δ⁺(A, i)) == sum(a::Arc -> t[a], δ⁻(A, i)) + sum(a::Arc -> x[a] * time(data, a), δ⁺(A, i)))
            @constraint(model, ub[i::Int in V′], sum(a::Arc -> t[a], δ⁻(A, i)) <= T - sum(block::Vi -> y[block] * blockTime(data, block), B))
            @constraint(model, ub1[i::Int in V′], sum(a::Arc -> t[a], δ⁻(A, i)) <= sum(a::Arc -> x[a], δ⁻(A, i)) * T)
            @constraint(model, ub2[a::Arc in A], t[a] <= x[a] * T)
        else # Nodes MTZ
            @variable(model, t[i::Int in keys(V)], lower_bound = 0, upper_bound = T)
            @constraint(model, t[depot] == 0.0)
            @constraint(model, mtz[a::Arc in A′], t[last(a)] >= t[first(a)] + x[a] * time(data, a) - (1 - x[a]) * T - x[reverse(a)] * time(data, reverse(a)))
            @constraint(model, ub1[i::Int in V′], t[i] <= T - sum(block::Vi -> y[block] * blockTime(data, block), B))
            @constraint(model, ub2[i::Int in V′], t[i] <= sum(a::Arc -> x[a], δ⁺(A, i)) * T)
            @constraint(model, ub3[i::Int in V′], t[i] <= sum(a::Arc -> x[a], δ⁻(A, i)) * T)
        end
        # subtour cuts
        addSubtourCuts(model, subtour_cuts) 
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
        model[:x], model[:y] = x, y

        return model
    end

    # creating model with both variables (x and y) relaxed
    model::Model = createModel(true, true)

    # getting intersection cuts
    if app["intersection-cuts"]
        @debug "Getting intersection cuts"
        
        elapsed_time::Float64 = @elapsed intersection_cuts1::Vector{Arcs}, intersection_cuts2::Vector{Arcs} = getIntersectionCuts(data) # get intersection cuts

        info["intersectionCutsTime"] = string(elapsed_time)
        info["intersectionCuts1"] = string(length(intersection_cuts1))
        info["intersectionCuts2"] = string(length(intersection_cuts2))

        # adding intersection cuts to the recently create model
        addIntersectionCuts1(model, intersection_cuts1)
        addIntersectionCuts2(model, intersection_cuts2)
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
    if app["subcycle-separation"]
        model = createModel(true, !app["y-integer"])
        
        optimize!(model)

        info["maxFlowCutsTime"] = string(@elapsed subtour_cuts = getSubtourCuts(data, model, info))
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
end
