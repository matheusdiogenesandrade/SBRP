
include("dd.jl")
include("dd_exact_test.jl")
include("dd_relaxed_test.jl")

using .DecisionDiagram
using .ExactDecisionDiagramTests
using .RelaxedDecisionDiagramTests

function main()
  RelaxedDecisionDiagramTests.test()
end
