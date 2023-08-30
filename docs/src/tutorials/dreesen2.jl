# # Getting rid of root at infinity
#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/dreesen2.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/dreesen2.ipynb)
# **Adapted from**: Section~2.2.1 p. 33 of [D13]
#
# [D13] Dreesen, Philippe.
# *Back to the Roots: Polynomial System Solving Using Linear Algebra*
# Ph.D. thesis (2013)

using LinearAlgebra
using TypedPolynomials
using Macaulay
using JuMP
using MultivariateMoments

# Consider the system given in [D13, (2.3)] which corresponds to `Systems/dreesen2` of [macaulaylab](http://www.macaulaylab.net/):

@polyvar x y z
system = [
    x^2 - x*y + z,
    2y^3 - 2x*y^2 - 3x * y,
    z^3 - x*y*z - 2,
]

# We first try to solve the system:

sols = solve_system(system, column_maxdegree = 6)
nothing # hide

# The real solutions are

sols

# ## Staicase analysis

solver = Iterator(system, Macaulay.Solver())
step!(solver)

# We can look at the solver as follows:

solver

# We can see in the border dependence that even if `x^2` is a corner, the multiplication
# matrix for `x` cannot be computed as `y^3` (resp. `x*z^2`, `y*z^2`) is standard but
# `x * y^3` (resp. `x^2 * z^2`, `x * y * z^2`) is indepedent.

using Plots
plot(saturated_dependence(solver))

# Let's do another step:

step!(solver)

# This time, the blocker for computing the multiplication matrix for `x`
# are `z^4` and `y * z^3` among others.

plot(saturated_dependence(solver))

# Let's do another step:

step!(solver)

# Now they are `z^5` and `y * z^4` among others.

plot(saturated_dependence(solver))

# Let's do a last step:

step!(solver)

# Now we see that the whole border is dependent so the four
# multiplication matrices can be computed.

plot(saturated_dependence(solver))

# In retrospect, we we probably should have expanded in priority towards larger
# exponents for `z`.

# ### Saturation
#
# Let's start again by starting with the same first stop.

solver = Iterator(system, Macaulay.Solver())
step!(solver)

# This time, let's focus on saturating the first
# non-saturated standard monomials:

step!(solver, FirstStandardNonSaturated(10))

# We can see that `z` is now saturated:

plot(saturated_dependence(solver))

# Let's saturate the next one:

step!(solver, FirstStandardNonSaturated(10))

# We can see that `y` went from trivial to saturated:

plot(saturated_dependence(solver))

# Let's saturate the next one:

step!(solver, FirstStandardNonSaturated(10))

# We can see that `x` went from trivial to saturated:

plot(saturated_dependence(solver))

# Let's saturate the next one:

step!(solver, FirstStandardNonSaturated(10))

# We can see that `z^2` is now saturated:

plot(saturated_dependence(solver))

# We can see that this does not seem to work.
# By saturating a monomial,
# we generate all the rows that have a nonzero entry for the corresponding
# column in the Macaulay matrix.
# If the monomial is standard, it means that the column is a linear combination of other columns
# (which are dependent in the Macaulay Nullspace)
# of the Macaulay matrix.
# If we saturate these other columns, the linear combination may not work anymore.
# But we are only saturating the standard monomials, we don't saturate
# the rows that are dependent in the Macaulay Nullspace which explains
# the result we have seen here.

# ## Moment approach

# With moment matrix of degree 3:

import SCS
solver = SCS.Optimizer
M = moment_matrix(system, solver, 3)
nothing # hide

atomic_measure(M, 1e-4, ShiftNullspace())

# With moment matrix of degree 4:

M = moment_matrix(system, solver, 4)
nothing # hide

atomic_measure(M, 1e-4, ShiftNullspace())

# With moment matrix of degree 5:

M = moment_matrix(system, solver, 5)
nothing # hide

atomic_measure(M, 1e-4, ShiftNullspace())
