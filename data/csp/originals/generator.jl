using Random

const SEED = 12345
const MAX_PROFIT = 70
const OUTPUT_DIR = "../"

Random.seed!(SEED);

instances = [
             "ali535-3.csp2",
             "ali535-5.csp2",
             "ali535-7.csp2",
             "att532-3.csp2",
             "att532-5.csp2",
             "att532-7.csp2",
             "berlin52-11.csp2",
             "berlin52-7.csp2",
             "berlin52-9.csp2",
             "d493-11.csp2",
             "d493-7.csp2",
             "d493-9.csp2",
             "d657-3.csp2",
             "d657-5.csp2",
             "d657-7.csp2",
             "eil51-11.csp2",
             "eil51-7.csp2",
             "eil51-9.csp2",
             "eil76-11.csp2",
             "eil76-7.csp2",
             "eil76-9.csp2",
             "fl417-11.csp2",
             "fl417-7.csp2",
             "fl417-9.csp2",
             "gil262-11.csp2",
             "gil262-7.csp2",
             "gil262-9.csp2",
             "gr666-3.csp2",
             "gr666-5.csp2",
             "gr666-7.csp2",
             "kroA100-11.csp2",
             "kroA100-7.csp2",
             "kroA100-9.csp2",
             "kroA150-11.csp2",
             "kroA150-7.csp2",
             "kroA150-9.csp2",
             "kroA200-11.csp2",
             "kroA200-7.csp2",
             "kroA200-9.csp2",
             "kroB100-11.csp2",
             "kroB100-7.csp2",
             "kroB100-9.csp2",
             "kroB150-11.csp2",
             "kroB150-7.csp2",
             "kroB150-9.csp2",
             "kroB200-11.csp2",
             "kroB200-7.csp2",
             "kroB200-9.csp2",
             "kroC100-11.csp2",
             "kroC100-7.csp2",
             "kroC100-9.csp2",
             "kroD100-11.csp2",
             "kroD100-7.csp2",
             "kroD100-9.csp2",
             "kroE100-11.csp2",
             "kroE100-7.csp2",
             "kroE100-9.csp2",
             "lin318-11.csp2",
             "lin318-7.csp2",
             "lin318-9.csp2",
             "p654-3.csp2",
             "p654-5.csp2",
             "p654-7.csp2",
             "pcb442-11.csp2",
             "pcb442-7.csp2",
             "pcb442-9.csp2",
             "pr226-11.csp2",
             "pr226-7.csp2",
             "pr226-9.csp2",
             "pr264-11.csp2",
             "pr264-7.csp2",
             "pr264-9.csp2",
             "pr299-11.csp2",
             "pr299-7.csp2",
             "pr299-9.csp2",
             "pr439-11.csp2",
             "pr439-7.csp2",
             "pr439-9.csp2",
             "pr76-11.csp2",
             "pr76-7.csp2",
             "pr76-9.csp2",
             "rat575-3.csp2",
             "rat575-5.csp2",
             "rat575-7.csp2",
             "rat783-3.csp2",
             "rat783-5.csp2",
             "rat783-7.csp2",
             "rat99-11.csp2",
             "rat99-7.csp2",
             "rat99-9.csp2",
             "rd100-11.csp2",
             "rd100-7.csp2",
             "rd100-9.csp2",
             "rd400-11.csp2",
             "rd400-7.csp2",
             "rd400-9.csp2",
             "st70-11.csp2",
             "st70-7.csp2",
             "st70-9.csp2",
             "ts225-11.csp2",
             "ts225-7.csp2",
             "ts225-9.csp2",
             "tsp225-11.csp2",
             "tsp225-7.csp2",
             "tsp225-9.csp2",
             "u574-3.csp2",
             "u574-5.csp2",
             "u574-7.csp2",
             "u724-3.csp2",
             "u724-5.csp2",
             "u724-7.csp2",
            ]

for instance in instances

    open(instance) do f::IOStream

        new_instance = ""
        # ignore the first line
        new_instance *= readline(f) * "\n"

        # get number of nodes
        line = readline(f)
        new_instance *= line * "\n"

        nNodes::Int = parse(Int, split(line, [' ']; limit=0, keepempty=false)[end])

        # ignore number of nodes per cluster
#        new_instance *= readline(f) * "\n"
        readline(f)

        # ignore distance matrix title
        new_instance *= readline(f) * "\n"

        # get distances
        for i in 1:nNodes
            new_instance *= readline(f) * "\n"
        end

        # get blocks

        # ignore clusters list title
        new_instance *= readline(f) * "\n"

        #
        clusters::Vector{Vector{Int}} = Vector{Vector{Int}}()
        for _ in 1:nNodes

            cluster::Vector{Int} = map(k::SubString{String} -> parse(Int, k) + 1, split(readline(f), [',', ' ']; limit=0, keepempty=false))
            push!(clusters, cluster)

        end

        blocks_set::Set{Vector{Int}} = Set{Vector{Int}}()
        for i::Int in 1:nNodes 

            block::Vector{Int} = sort(filter(j::Int -> i in clusters[j], 1:nNodes))

            push!(blocks_set, sort(block))

        end

        # define profit
        for block::Vector{Int} in blocks_set
            new_instance *= join(block, " ") * " " * string(rand(1:MAX_PROFIT)) * "\n"
        end

        new_instance *= readline(f) * "\n"

        open(OUTPUT_DIR * instance, "w") do file
            write(file, new_instance)
        end
    end

end
