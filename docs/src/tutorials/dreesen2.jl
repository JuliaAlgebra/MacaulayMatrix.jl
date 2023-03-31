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

sols = solve_system(system, 6)
nothing # hide

# The real solutions are

sols

# With moment matrix of degree 3:

import SCS
solver = SCS.Optimizer
M = moment_matrix(system, solver, 3)
extractatoms(M, 1e-4, ShiftNullspace())

# With moment matrix of degree 4:

M = moment_matrix(system, solver, 4)
extractatoms(M, 1e-4, ShiftNullspace())

# With moment matrix of degree 5:

M = moment_matrix(system, solver, 5)
extractatoms(M, 1e-4, ShiftNullspace())
