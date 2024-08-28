using LinearAlgebra
using MacaulayMatrix
using JuMP
using MultivariateMoments
using DynamicPolynomials
using HomotopyContinuation

# --------------- Example Benoit -------------------
DynamicPolynomials.@polyvar x[1:2]
sys = [x[2]^2 - x[1], x[2]^2]

# Real-valued solutions using Macaulay:
sols = solve_system(sys, column_maxdegree = 4, print_level = 3)
# All solutions using HomotopyContinuation:
HomotopyContinuation.solve(sys)

# ---------------- Example alternative ------------------
# [(x[1] - 1) * (x[2] - 3), (x[1] - 4) * (x[2] - 2im) * (x[2] + 2im)]
sys = [3 - x[2] - 3x[1] + x[1]x[2], -16 + 4x[1] - 4x[2]^2 + x[1]x[2]^2]
# Real-valued solutions using Macaulay:
sols = solve_system(sys, column_maxdegree = 4, print_level = 3)
# All solutions using HomotopyContinuation:
res = HomotopyContinuation.solve(sys)
results(res)

# The system has 3 isolated common roots (1, ±2im), (4, 3):
subs.(sys, x[1] => 1, x[2] => 2im)
subs.(sys, x[1] => 1, x[2] => -2im)
subs.(sys, x[1] => 4, x[2] => 3)

# Check whether nullity has stabilized for d = 4:
Z2 = nullspace(macaulay(sys, 2)).matrix
Z3 = nullspace(macaulay(sys, 3)).matrix
Z4 = nullspace(macaulay(sys, 4)).matrix
Z5 = nullspace(macaulay(sys, 5)).matrix
# Nullity stabilizes at 4. Bezout: 6?

# Check for degree gap:
~, S = svd(Z4[1:length(monomials(x, [0, 1])), :])
~, S = svd(Z4[1:length(monomials(x, [0, 1, 2])), :])
# Second degree block does not contain linearly independent row.
~, S = svd(Z4[1:length(monomials(x, [0, 1, 2, 3])), :])
~, S = svd(Z4)

# Moment approach:
d_max = 4
Z = nullspace(macaulay(sys, d_max))
d_m = 2
ν = moment_matrix(Z, big_clarabel, d_m, T = BigFloat)

# Manually inspect rank of Moment matrix via SVD:
M = value_matrix(ν)
M, S = svd(M)
round.(log10.(S))

# Given the rank r, check for solutions:
r = 3
res = atomic_measure(ν, FixedRank(r), ShiftNullspace())
