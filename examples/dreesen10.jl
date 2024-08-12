# # Positive dimensional roots at infinity
# **Adapted from**: Example 6.14 p. 102 of [D13]
#
# [D13] Dreesen, Philippe.
# *Back to the Roots: Polynomial System Solving Using Linear Algebra*
# Ph.D. thesis (2013)

using Test
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

expected(T=Float64) = [
    T[1//2, 1//2, √T(2//3), -√T(2//3)],
    T[1//2, 1//2, -√T(2//3), √T(2//3)],
]

# With the classical MacaulayMatrix approach, the nullity increases by 2 at every degree
# because of the positive dimensional solution set at infinity.
# The solution is obtained at degree `6`.

solve_system(system, column_maxdegree = 8)

# With moment matrix of degree 4:

include("solvers.jl")
big_clarabel = clarabel_optimizer(T = BigFloat)

ν4 = moment_matrix(system, big_clarabel, 4, T = BigFloat)

# The rank of the moment matrix is `9`:

r4 = MacaulayMatrix.cheat_rank(ν4, x, expected(BigFloat), LeadingRelativeRankTol(1e-6))

# The nullspace is then of degree 6 and we can see below that they are
# 6 valid equations for our solutions.

MacaulayMatrix.errors(ν4, x, expected(BigFloat))

# Using the nullspace from the SVD of the moment matrix, we get
# these 6 equations:

M4 = nullspace(ν4, FixedRank(r4))

# When we compare to the initial system, we see that these
# are simply redundant equations:

Δ4 = nonredundant(M4, system)

# We can also use the image space from the SVD. This is a
# Macaulay nullspace. We can therefore find there the dependent
# and independent monomials and put this nullspace in standard form.
# The nullspace of this standard form gives a Macaulay matrix
# that's also in standard form hence is a bit more readable.

M4 = nullspace(ν4, FixedRank(r4), ShiftNullspace())

# We lose one equation here, probably because we ignore dependent
# monomials that are multiple of another dependent one.
# Taking the difference, unsurprisingly gives us the same result,
# no new equation.

Δ4 = nonredundant(M4, system)

# With moment matrix of degree 4 from a Macaulay nullspace of degree 5,
# we get:

ν5 = moment_matrix(system, big_clarabel, 5, T = BigFloat)

# The rank of the moment matrix is now `5`:

r5 = MacaulayMatrix.cheat_rank(ν5, x, expected(BigFloat), LeadingRelativeRankTol(1e-6))

# Indeed, we see that the SVD gives 10 valid equations:

MacaulayMatrix.errors(ν5, x, expected(BigFloat))

# Among these, we have 6 new equations:

Δ5 = nonredundant(nullspace(ν5, FixedRank(r5), ShiftNullspace()), system)

# Are these equations also valid for complex solutions ?
