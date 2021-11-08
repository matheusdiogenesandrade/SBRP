using JuMP
using CPLEX

function solve(model)
  optimize!(model)
  # termination status MOI.OPTIMAL and 
  return termination_status(model) in [MOI.OPTIMAL, MOI.TIME_LIMIT]
end
