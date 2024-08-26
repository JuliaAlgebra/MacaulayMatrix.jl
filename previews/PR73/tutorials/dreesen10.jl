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
using MacaulayMatrix
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

# With the classical MacaulayMatrix approach, the nullity increases by 2 at every degree
# because of the positive dimensional solution set at infinity.

sols = #src
solve_system(system, column_maxdegree = 6)
@test length(sols) == 2 #src

# Let's compute the moment matrix of degree 4 from the nullspace of
# the Macaulay matrix of degree 5 using the Clarabel SDP solver:

import Clarabel
ν = moment_matrix(system, Clarabel.Optimizer, 5)

# The nullspace of this moment matrix is a Macaulay matrix

M = nullspace(ν, ShiftNullspace())

# Each coefficient is close to a rational number, after some clean up,
# we get the following system:

real_system = clean(M; tol = 1e-6).polynomials
@test in(2x[1] - 1, real_system) #src
@test in(2x[2] - 1, real_system) #src
@test in(x[3] + x[4], real_system) #src
@test in(2 + 3x[3] * x[4], real_system) #src

# We have the equations `2x[1] = 1` and `2x[2] = 1` from which we get
# `x[1] = 1/2` and `x[2] = 1/2`.
# We then have the equations `x[3] = -x[4]` and `3x[3] * x[4] = -2` from
# which we get `x[3] = ±√(2/3)` and `x[4] = ∓√(2/3)`.
# Notice that the complex solutions are not solutions of the system anymore.
# The Macaulay solver indeed find these real solutions direction at degree 2.

sols = #src
solve_system(real_system, column_maxdegree = 2)
@test length(sols) == 2 #src
