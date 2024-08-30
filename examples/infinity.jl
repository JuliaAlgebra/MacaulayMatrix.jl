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
sols = solve_system(
    sys,
    column_maxdegree = 4,
    sparse_columns = false,
    print_level = 3,
)
# All solutions using HomotopyContinuation:
HomotopyContinuation.solve(sys)

# ---------------- Example PhD Dreesen (Section 2.3.1) ------------------
DynamicPolynomials.@polyvar x[1:2]
sys = [3 − x[2] − 3x[1] + x[1]x[2], −16 + 4x[1] − 4x[2]^2 + x[1]x[2]^2]

# Real-valued solutions using Macaulay:
sols = solve_system(
    sys,
    column_maxdegree = 4,
    sparse_columns = false,
    print_level = 3,
)
# All solutions using HomotopyContinuation:
res = HomotopyContinuation.solve(sys)
results(res)

# 3 Affine solutions, 1 real-valued.
subs.(sys, x[1] => 4, x[2] => 3)
subs.(sys, x[1] => 1, x[2] => 2im)
subs.(sys, x[1] => 1, x[2] => -2im)

# Inspect nullity of Macaulay matrices:
Z2 = nullspace(macaulay(sys, 2, sparse_columns = false)).matrix
Z3 = nullspace(macaulay(sys, 3, sparse_columns = false)).matrix
Z4 = nullspace(macaulay(sys, 4, sparse_columns = false)).matrix
Z5 = nullspace(macaulay(sys, 5, sparse_columns = false)).matrix
Z6 = nullspace(macaulay(sys, 6, sparse_columns = false)).matrix
# Nullity stabilizes at 6. This corresponds to Bézout's thm. 
# There are 3 solutions at infinity. 

~, S = svd(Z4[1:length(monomials(x, [0, 1])), :])
~, S = svd(Z4[1:length(monomials(x, [0, 1, 2])), :])
~, S = svd(Z4[1:length(monomials(x, [0, 1, 2, 3])), :])
~, S = svd(Z4)
# second degree block does not contain linearly independent row(s).

# Is the Macaulay matrix 2-shift-invariant?
~, S = svd(Z4[1:length(monomials(x, [0, 1, 2])), :])
~, S = svd(Z5[1:length(monomials(x, [0, 1, 2])), :])
# Yes.

# Apply moment-matrix approach:
include("solvers.jl")
big_clarabel = clarabel_optimizer(T = BigFloat)

d_mac = 4
Z = nullspace(macaulay(sys, d_max, sparse_columns = false))
d_m = 1 # d_m = d_mom / 2
ν = moment_matrix(Z, big_clarabel, d_m, T = BigFloat)

# Manually inspect rank of Moment matrix via SVD:
M = value_matrix(ν)
M, S = svd(M)
round.(log10.(S))

# Given the rank r, check for solutions:
r = 1
res = atomic_measure(ν, FixedRank(r), ShiftNullspace())

# Check for flatness:
ν_trunc = truncate(ν, 0)
M = value_matrix(ν_trunc)
M, S = svd(M)
round.(log10.(S))

# ----------- What if we do not truncate roots at infinity? --------
d_mac = 4
Z = nullspace(macaulay(sys, d_max, sparse_columns = false))
d_m = 2 # d_m = d_mom / 2
ν = moment_matrix(Z, big_clarabel, d_m, T = BigFloat)

# Manually inspect rank of Moment matrix via SVD:
M = value_matrix(ν)
M, S = svd(M)
round.(log10.(S))

# Given the rank r, check for solutions:
r = 3
res = atomic_measure(ν, FixedRank(r), ShiftNullspace())
# Multiplication matrices have relatively large commutation error!
# So not possible to retrieve correct solutions. 

# Let's try truncating:
ν_t1 = truncate(ν, 1)
M = value_matrix(ν_t1)
M, S = svd(M)
round.(log10.(S))

# Given the rank r, check for solutions:
r = 1
res = atomic_measure(ν_t1, FixedRank(r), ShiftNullspace())
# Correct real-valued solution is identified! 

ν_t2 = truncate(ν, 0)
M = value_matrix(ν_t2)
M, S = svd(M)
round.(log10.(S))

# Given the rank r, check for solutions:
r = 1
res = atomic_measure(ν_t2, FixedRank(r), ShiftNullspace())