using TypedPolynomials
using SemialgebraicSets
import DynamicPolynomials
function system(k, solver)
    @polyvar x y
    return @set x == y^k && y^k == 0 solver
end
using Random
Random.seed!(80)

using MacaulayMatrix
sys_gap = system(4, MacaulayMatrix.Solver(sparse_columns = false, wait_for_gap = true))
solutions = collect(sys_gap)
nothing #hide

solutions

using MacaulayMatrix
sys = system(4, MacaulayMatrix.Solver(sparse_columns = false))
solutions = collect(sys)
nothing #hide

solutions

sys_sparse = system(4, MacaulayMatrix.Solver())
solutions = collect(sys_sparse)
nothing #hide

solutions

solver_gap = init(sys_gap, sys_gap.solver)

step!(solver_gap)

step!(solver_gap)

step!(solver_gap)

step!(solver_gap)

solver_gap

using Plots
plot(saturated_dependence(solver_gap))

solver_gap = init(sys_gap, MacaulayMatrix.Solver(dependence = MacaulayMatrix.MM.LinearDependence, sparse_columns = false, wait_for_gap = true))

solve!(solver_gap)

using Plots
plot(saturated_dependence(solver_gap))

solver = init(sys, sys.solver)

step!(solver)

solver

plot(saturated_dependence(solver))

sparse_solver = init(sys_sparse, sys_sparse.solver)

step!(sparse_solver)

sparse_solver

plot(saturated_dependence(sparse_solver))

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
