# # Positive dimensional roots at infinity
#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/dreesen10.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/dreesen10.ipynb)
# **Adapted from**: Example 6.14 p. 102 of [D13]
#
# [D13] Dreesen, Philippe.
# *Back to the Roots: Polynomial System Solving Using Linear Algebra*
# Ph.D. thesis (2013)

using Test #src
using LinearAlgebra
using TypedPolynomials
using Macaulay
using JuMP
using MultivariateMoments

# Consider the system given in [D13, Example 6.14] which corresponds to `Systems/dreesen10` of [macaulaylab](http://www.macaulaylab.net/):

@polyvar x[1:4]
system = [
    x[1] + x[2] - 1,
    x[1] * x[3] + x[2] * x[4],
    x[1] * x[3]^2 + x[2] * x[4]^2 - 2/3,
    x[1] * x[3]^3 + x[2] * x[4]^3,
]

# With the classical Macaulay approach, the nullity increases by 2 at every degree
# because of the positive dimensional solution set at infinity.

sols = solve_system(system, 8)

# With moment matrix of degree 6:

import SCS
solver = SCS.Optimizer
M = moment_matrix(system, solver, 6)

# We don't find anything:

extractatoms(M, 1e-4, ShiftNullspace())

# With moment matrix of degree 7:

M = moment_matrix(system, solver, 7)

# We get different solutions for different runs because of the random combinations of multiplication matrices:
# so let's fix the seed to make it reproducible:

using Random
Random.seed!(0)
extractatoms(M, 1e-4, ShiftNullspace())

# The second time, no luck:

extractatoms(M, 1e-4, ShiftNullspace())

# The third one contains the solutions `(0.5, 0.5, 0.8165, -0.8165)`
# and `(0.5, 0.5, -0.8165, 0.8165)`:

sols = extractatoms(M, 1e-4, ShiftNullspace())
@test check(sols, [0.5, 0.5, -0.81, 0.81]) #src
@test check(sols, [0.5, 0.5, 0.81, -0.81]) #src
