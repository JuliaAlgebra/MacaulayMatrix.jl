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
compute_support!(ν4, FixedRank(r4), ImageSpaceSolver(SVDLDLT(), Echelon()))
support_error(ν4, x, expected(BigFloat))
M4 = nullspace(ν4)
M4 = nullspace(ν4, FixedRank(r4))
solve_system([system; Float64.(M4) * ν4.basis.monomials])

ν5 = moment_matrix(system, big_clarabel, 5, T = BigFloat)
cheat_system(ν5, x, expected(BigFloat), LeadingRelativeRankTol(1e-6))
MacaulayMatrix.errors(ν5, x, expected(BigFloat))
r5 = MacaulayMatrix.cheat_rank(ν5, x, expected(BigFloat), LeadingRelativeRankTol(1e-6))
compute_support!(ν5, FixedRank(r5), ImageSpaceSolver(SVDLDLT(), Echelon()))
ν5.support
support_error(ν5, x, expected(BigFloat))
M5 = nullspace(ν5)
M5 * ν5.basis.monomials
solve_system(M5 * ν5.basis.monomials)
solve_system(Float64.(M5) * ν5.basis.monomials)

atomic_measure(ν5, FixedRank(r5))
atomic_measure(ν5, FixedRank(r5), ShiftNullspace())

ν6 = moment_matrix(system, solver, 6, T = BigFloat)
solutions(ν6)
compute_support!(ν6, FixedRank(25))
ν6.support
using SemialgebraicSets
for eq in equalities(ν6.support)
    @show eq(x => [0.5, 0.5, -0.8165, 0.8165])
end

# We don't find anything:

atomic_measure(M, 1e-5, ShiftNullspace())

# With moment matrix of degree 7:

ν7 = moment_matrix(system, solver, 7)
solutions(ν7)

ν8 = moment_matrix(system, solver, 8)
solutions(ν8)

# We get different solutions for different runs because of the random combinations of multiplication matrices:
# so let's fix the seed to make it reproducible:

using Random
Random.seed!(0)
atomic_measure(M, 1e-4, ShiftNullspace())

# The second time, no luck:

atomic_measure(M, 1e-4, ShiftNullspace())

# The third one contains the solutions `(0.5, 0.5, 0.8165, -0.8165)`
# and `(0.5, 0.5, -0.8165, 0.8165)`:

sols = atomic_measure(M, 1e-6, ShiftNullspace())
check(sols, x) = any(atom -> isapprox(atom.center, x, rtol=1e-2), sols.atoms)
@test check(sols, [0.5, 0.5, -0.81, 0.81])
@test check(sols, [0.5, 0.5, 0.81, -0.81])

display(solutions(M))
c = sols[1]

M = moment_matrix(system, solver, 8)
display(solutions(M))
