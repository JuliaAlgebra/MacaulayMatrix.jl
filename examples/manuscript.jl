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

# TODO: solve_system() removes complex-valued roots! How to compute all of them?
sols = solve_system(sys, column_maxdegree = 4, print_level = 3)

# Let's try the moment-matrix approach:
d_max = 4
Z = nullspace(macaulay(sys, d_max))

include("solvers.jl")
big_clarabel = clarabel_optimizer(T = BigFloat)

M = moment_matrix(Z, big_clarabel, 2, T = BigFloat)
res = atomic_measure(M, 1e-4, ShiftNullspace())
# The solution seems off. We probably need to be more carefull with the rank checks:

# Let's cheat to find the correct rank:
expected(T = Float64) = [T[1, 4], T[3, 3]]
r = MacaulayMatrix.cheat_rank(
    M,
    [x, y],
    expected(BigFloat),
    LeadingRelativeRankTol(1e-6),
)

# Now try again:
res = atomic_measure(M, FixedRank(r), ShiftNullspace())