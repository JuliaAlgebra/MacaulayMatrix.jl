# # Degree gap versus staircase gap
#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/infinity.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/infinity.ipynb)

using TypedPolynomials
using SemialgebraicSets
import DynamicPolynomials
function system(k, solver)
    @polyvar x y
    return @set x == y^k && y^k == 0 solver
end
using Random
Random.seed!(80)

# Let's start with the classical approach with all columns up to degree 5 and waiting for the gap:

using Macaulay
sys_gap = system(4, Macaulay.Solver(default_iteration = ColumnDegreeIteration(sparse_columns = false, wait_for_gap = true)))
collect(sys_gap)

# If we don't wait for the gap, we get it earlier:

using Macaulay
sys = system(4, Macaulay.Solver(default_iteration = ColumnDegreeIteration(sparse_columns = false)))
collect(sys)

# But we can actually also just use the 2 columns that are actually used:

sys_sparse = system(4, Macaulay.Solver(default_iteration = ColumnDegreeIteration()))
collect(sys_sparse)

# ## In more details

solver_gap = init(sys_gap, sys_gap.solver)

# After one step, we don't have the solution yet

step!(solver_gap)

# After the second step, still no solution

step!(solver_gap)

# After the third step, still no solution

step!(solver_gap)

# After the fourth step, we find the solution

step!(solver_gap)

# We can inspect the solver at follows

solver_gap

# The border dependence can be plotted for more visual inspection:

using Plots
plot(solver_gap.border.dependence)

# We can see that the monomials of larger degree have not been used here.
# We can see them used with `AnyDependence` instead of `StaircaseDependence`

solver_gap = init(sys_gap, Macaulay.Solver(dependence = Macaulay.MM.AnyDependence, default_iteration = ColumnDegreeIteration(sparse_columns = false, wait_for_gap = true)))

# Let's solve again but not all the step at once:

solve!(solver_gap)

# This time, we see the monomials at infinity as the blue balls outside the staircase.
# We can see that the degree 4 contains no independent as it is the gap zone.

using Plots
plot(solver_gap.border.dependence)

# ## Not waiting for the gap in more details

solver = init(sys, sys.solver)

# After one step, we find the solution

step!(solver)

# We can inspect `m` for more details

solver

# The border dependence can be plotted for more visual inspection:

plot(solver.border.dependence)

# Even if there is no gap, the border is complete so we can get the multiplication matrices.

# ## Sparse columns in more details

sparse_solver = init(sys_sparse, sys_sparse.solver)

# After one step, we find the solution

step!(sparse_solver)

# We can inspect `m` for more details

sparse_solver

# The border dependence can be plotted for more visual inspection:

plot(sparse_solver.border.dependence)

# As we can see, for the sparse one, the standard monomials are "trivial"
# because they are trivially detected as independent since they are not part of the basis.
# The border is missing but the multiplication matrix for `y` can be computed first
# and then, using this multiplication matrix, the missing relations for the border
# can be obtained and then the multiplication matrix for `x` can be obtained.
