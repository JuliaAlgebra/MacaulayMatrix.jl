# ---------- All numerical examples used throughout the manuscript --------------

using LinearAlgebra
using MacaulayMatrix
using JuMP
using MultivariateMoments
using DynamicPolynomials

# Removing complex-valued roots:
DynamicPolynomials.@polyvar x y
sys = [4 - 5x + x^2, -12 + 4 * y - 3 * y^2 + y^3]

# The system has 6 isolated common roots (1, ±2im), (4, ±2im) and (1, 3), (4, 3):
subs.(sys, x => 1, y => 2im)
subs.(sys, x => 1, y => -2im)
subs.(sys, x => 4, y => 2im)
subs.(sys, x => 4, y => -2im)
subs.(sys, x => 1, y => 3)
subs.(sys, x => 4, y => 3)

# ----------- Standard Macaulay nullspace approach ---------------

# TODO: solve_system() removes complex-valued roots! How to compute all of them?
sols = solve_system(sys, column_maxdegree = 4, print_level = 3)

# Check whether nullity has stabilized for d = 4:
Z3 = nullspace(macaulay(sys, 3)).matrix
Z4 = nullspace(macaulay(sys, 4)).matrix
Z5 = nullspace(macaulay(sys, 5)).matrix

# Nullity stabilizes from d = 4, ...
# Check for degree gap:
~, S = svd(Z4[1:length(monomials((x, y), [0, 1, 2])), :])
~, S = svd(Z4[1:length(monomials((x, y), [0, 1, 2, 3])), :])
~, S = svd(Z4)
# There is a degree gap: the final degree block (d = 4) does not contain linearly independent rows.

# ----------- Automated moment approach ---------------

# Let's try the moment-matrix approach:
d_max = 4
Z = nullspace(macaulay(sys, d_max))

include("solvers.jl")
big_clarabel = clarabel_optimizer(T = BigFloat)

ν = moment_matrix(Z, big_clarabel, 2, T = BigFloat)
res = atomic_measure(ν, 1e-4, ShiftNullspace())
# We are able to retrieve the two real-valued solutions!

# ------------ Manual hyperparameter search -----------
# Hyperparameters:
# 1) d_max = d_mac: degree of Macaulay to compute nullspace Z
# 2) d_m = d_mom: degree of the PSD Moment matrix M to construct from Z (d_m <= d_max)
# 3) r: rank of the Moment matrix M (used to retrieve the real-valued solutions)

# Let's try the moment-matrix approach:
d_max = 4
Z = nullspace(macaulay(sys, d_max))
d_m = 2
ν = moment_matrix(Z, big_clarabel, d_m, T = BigFloat)

# Manually inspect rank of Moment matrix via SVD:
M = value_matrix(ν)
M, S = svd(M)
round.(log10.(S))

# Given the rank r, check for solutions:
r = 2
res = atomic_measure(ν, FixedRank(r), ShiftNullspace())

# Results (Clarabel BigFloat):
# d_max = 4, d_m = 1: rank is not clear, seems full rank (3)
# d_max = 4, d_m = 2: rank is clearly 2, two valid solutions! 
# d_max = 4, d_m = 3: rank is clearly 6, no solutions... 

# d_max = 5, d_m = 1: full rank
# d_max = 5, d_m = 2: rank is clearly 2, two valid solutions!
# d_max = 5, d_m = 3: rank is clearly 6, two valid solutions! 
# d_max = 5, d_m = 4: rank seems 7, not possible to find solutions.
# d_max = 5, d_m = 5: difficult to find rank gap. 

# d_max = 6, d_m = 1: full rank
# d_max = 6, d_m = 2: rank is clearly 2, two valid solutions!
# d_max = 6, d_m = 3: rank is clearly 2, two valid solutions!
# d_max = 6, d_m = 4: rank seems 7, two valid solutions! 

# d_max = 6, d_m = 5: rank = 7 or rank = 8, two valid solutions!

# TODO: when do we have `flatness'?

# -------------- `CheatRank' function --------------

d_max = 4
Z = nullspace(macaulay(sys, d_max))
d_m = 2
ν = moment_matrix(Z, big_clarabel, d_m, T = BigFloat)

# Let's cheat to find the correct rank:
expected(T = Float64) = [T[1, 3], T[4, 3]]
r = MacaulayMatrix.cheat_rank(
    ν,
    [x, y],
    expected(BigFloat),
    LeadingRelativeRankTol(1e-6),
)
res = atomic_measure(ν, FixedRank(r), ShiftNullspace())

# --------------- Truncate the obtained moment matrix ---------

d_max = 4
Z = nullspace(macaulay(sys, d_max))
d_m = 2
ν = moment_matrix(Z, big_clarabel, d_m, T = BigFloat)
ν_1 = truncate(ν, 1)

M = value_matrix(ν_1)
M, S = svd(M)
round.(log10.(S))

# Given the rank r, check for solutions:
r = 2
res = atomic_measure(ν_1, FixedRank(r), ShiftNullspace())
# Two real-valued solutions! 

# Alternative:
d_m = 1
ν = moment_matrix(Z, big_clarabel, d_m, T = BigFloat)
M = value_matrix(ν)
M, S = svd(M)
round.(log10.(S))

# Given the rank r, check for solutions:
r = 2
res = atomic_measure(ν, FixedRank(r), ShiftNullspace())

# --------------- Real radical ideal ---------------

M = nullspace(ν, ShiftNullspace(), FixedRank(r))
Δ = nonredundant(M, sys)
# TODO: get this to work.