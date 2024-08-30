# ---------- All numerical examples used throughout the manuscript --------------

using LinearAlgebra
using MacaulayMatrix
using JuMP
using MultivariateMoments
using DynamicPolynomials

# ---------------- Example ------------------
DynamicPolynomials.@polyvar x[1:2]
sys = [4 - 5x[1] + x[1]^2, -12 + 4 * x[2] - 3 * x[2]^2 + x[2]^3]

# The system has 6 isolated common roots (1, ±2im), (4, ±2im) and (1, 3), (4, 3):
subs.(sys, x[1] => 1, x[2] => 2im)
subs.(sys, x[1] => 1, x[2] => -2im)
subs.(sys, x[1] => 4, x[2] => 2im)
subs.(sys, x[1] => 4, x[2] => -2im)
subs.(sys, x[1] => 1, x[2] => 3)
subs.(sys, x[1] => 4, x[2] => 3)

# Real-valued solutions using Macaulay:
sols = solve_system(
    sys,
    column_maxdegree = 5,
    sparse_columns = false,
    print_level = 3,
)
# All solutions using HomotopyContinuation:
res = HomotopyContinuation.solve(sys)
results(res)

# Check whether nullity has stabilized for d = 4:
Z2 = nullspace(macaulay(sys, 2, sparse_columns = false)).matrix
Z3 = nullspace(macaulay(sys, 3, sparse_columns = false)).matrix
Z4 = nullspace(macaulay(sys, 4, sparse_columns = false)).matrix
Z5 = nullspace(macaulay(sys, 5, sparse_columns = false)).matrix
# Nullity stabilizes at 6. This corresponds to Bézout's thm.

# Check for degree gap:
~, S = svd(Z3[1:length(monomials(x, [0, 1])), :])
~, S = svd(Z3[1:length(monomials(x, [0, 1, 2])), :])
~, S = svd(Z3)
# No degree gap... 

~, S = svd(Z4[1:length(monomials(x, [0, 1])), :])
~, S = svd(Z4[1:length(monomials(x, [0, 1, 2])), :])
~, S = svd(Z4[1:length(monomials(x, [0, 1, 2, 3])), :])
~, S = svd(Z4)
# Fourth degree block does not contain linearly independent row.

# Is the Macaulay matrix 3-shift-invariant?
~, S = svd(Z4[1:length(monomials(x, [0, 1, 2, 3, 4])), :])
~, S = svd(Z5[1:length(monomials(x, [0, 1, 2, 3, 4])), :])
# Yes.

# ------------------ Moment approach ----------------
include("solvers.jl")
big_clarabel = clarabel_optimizer(T = BigFloat)

d_mac = 4
Z = nullspace(macaulay(sys, d_mac, sparse_columns = false))
d_m = 2 # d_m = d_mom / 2
ν = moment_matrix(Z, big_clarabel, d_m, T = BigFloat)

# Let's try the automated rank-check:
res = atomic_measure(ν, 1e-4, ShiftNullspace())
# We are able to retrieve the two real-valued solutions!

# Manually inspect rank of Moment matrix via SVD:
M = value_matrix(ν)
M, S = svd(M)
round.(log10.(S))

# Given the rank r, check for solutions:
r = 2
res = atomic_measure(ν, FixedRank(r), ShiftNullspace())
# Attention: normally you should truncate to remove effects of roots at infinity!

# Results (Clarabel BigFloat):
# d_mac = 4, d_mom = 1: rank is not clear, seems full rank (3)
# d_mac = 4, d_mom = 2: rank is clearly 2, two valid solutions! 

# d_mac = 5, d_mom = 1: full rank
# d_mac = 5, d_mom = 2: rank is clearly 2, two valid solutions!

# d_mac = 6, d_mom = 1: full rank
# d_mac = 6, d_mom = 2: rank is clearly 2, two valid solutions!
# d_mac = 6, d_mom = 3: rank is clearly 2, two valid solutions!

# --------------- Truncate the obtained moment matrix ---------

# `flatness' -> rank of two truncations does not increase

# Lets try (d_mac, d_mom, d_tc) = (4, 4, 2):
d_mac = 4
Z = nullspace(macaulay(sys, d_mac))
d_m = 2
ν = moment_matrix(Z, big_clarabel, d_m, T = BigFloat)
ν_trunc = truncate(ν, 1)

M = value_matrix(ν_trunc)
M, S = svd(M)
round.(log10.(S))

# Given the rank r, check for solutions:
r = 2
res = atomic_measure(ν_trunc, FixedRank(r), ShiftNullspace())

# Alternative, lets try (d_mac, d_mom, d_tc) = (4, 2, 2):
d_m = 1
ν_1 = moment_matrix(Z, big_clarabel, d_m, T = BigFloat)
M = value_matrix(ν_1)
M, S = svd(M)
round.(log10.(S))

# Given the rank r, check for solutions:
r = 2
res = atomic_measure(ν_1, FixedRank(r), ShiftNullspace())

# -------------- `CheatRank' function --------------

d_mac = 4
Z = nullspace(macaulay(sys, d_mac))
d_m = 2
ν = moment_matrix(Z, big_clarabel, d_m, T = BigFloat)

# Let's cheat to find the correct rank:
expected(T = Float64) = [T[4, 3]]
r = MacaulayMatrix.cheat_rank(
    ν,
    [x, y],
    expected(BigFloat),
    LeadingRelativeRankTol(1e-6),
)
res = atomic_measure(ν, FixedRank(r), ShiftNullspace())