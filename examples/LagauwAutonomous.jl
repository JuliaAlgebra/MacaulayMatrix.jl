using LinearAlgebra
using TypedPolynomials
using MacaulayMatrix
using JuMP
using MultivariateMoments
using DynamicPolynomials

# Walsh MEP from (Lagauw, Vanpoucke, De Moor, 2024: "Exact characterization of the global optima of least squares realization of autonomous LTI models as a multiparameter eigenvalue problem")
function composeSystemWalsh(y::Vector{Float64}, n::Integer)
    DynamicPolynomials.@polyvar a[1:n]
    aa = [1; a]

    N = length(y)
    T = zeros(eltype(aa), N - n, N)
    for i in 1:N-n
        T[i, i:i-1+length(aa)] = reverse(aa)
    end
    S = zeros(eltype(aa), N - 2n, N - n)
    for i in 1:N-2n
        S[i, i:i-1+length(aa)] = reverse(aa)
    end

    DynamicPolynomials.@polyvar g[1:N-2n]

    system = Vector{polynomial_type(eltype(aa), Float64)}(undef, N - n)
    system[1:N-n] = T * y - T * transpose(T) * transpose(S) * g

    return system
end

n = 1;
y = [4, 3, 2, 1];
y = convert(Vector{Float64}, y)
# D = 3 (degree of polynomials is at most cubic)
sys = composeSystemWalsh(y, n)

# TypedPolynomials.@polyvar a g[1:2]
# system = [
#     3 + 4 * a - (a * (g[1] + g[2] * a) + g[1] * a + g[1] * a^3)
#     2 + 3 * a - (g[1] + a * (a * (g[1] + g[2] * a) + g[1] * a) + 2 * g[2] * a)
#     1 + 2 * a - (g[2] + a * (g[1] + 2 * g[2] * a))
# ]

# ---------- Standard Macaulay Null space ------------

# Degree 7 Macaulay matrix suffices to retrieve V_r(I), by selection from V_c(I):
sols = solve_system(sys, column_maxdegree = 8, print_level = 3)

# Compute all solutions:
using HomotopyContinuation
res = HomotopyContinuation.solve(sys)
results(res)
# 6 complex-valued solutions, 1 real-valued solution.

# Step by step analysis:
solver = Iterator(sys, MacaulayMatrix.Solver())

# first shift, d: D --> D+1   
step!(solver)
solver
# d: D+1 --> D+2   
step!(solver)
solver
# d: D+2 --> D+3   
step!(solver)
solver
# d: D+3 --> D+4   
step!(solver)
solver

# Inspect shift-invariance of Macaulay:
Z3 = nullspace(macaulay(sys, 3, sparse_columns = false)).matrix
Z4 = nullspace(macaulay(sys, 4, sparse_columns = false)).matrix
Z5 = nullspace(macaulay(sys, 5, sparse_columns = false)).matrix
Z6 = nullspace(macaulay(sys, 6, sparse_columns = false)).matrix
Z7 = nullspace(macaulay(sys, 7, sparse_columns = false)).matrix
Z8 = nullspace(macaulay(sys, 8, sparse_columns = false)).matrix
Z9 = nullspace(macaulay(sys, 9, sparse_columns = false)).matrix
# Nullity stabilizes has constant rank increases of 5: solution set is 1-dimensional.

# Check for degree gap:
~, S = svd(Z7[1:length(monomials(x, [0, 1])), :])
# r = 5
~, S = svd(Z7[1:length(monomials(x, [0, 1, 2])), :])
# r = 7
~, S = svd(Z7[1:length(monomials(x, [0, 1, 2, 3])), :])
# r = 11
~, S = svd(Z7[1:length(monomials(x, [0, 1, 2, 3, 4])), :])
# r > 11
# ... 
~, S = svd(Z7)
# No degree gap.

~, S = svd(Z8[1:length(monomials(x, [0, 1])), :])
# r = 5
~, S = svd(Z8[1:length(monomials(x, [0, 1, 2])), :])
# r = 7
~, S = svd(Z8[1:length(monomials(x, [0, 1, 2, 3])), :])
# r = 7 --> GAP!
~, S = svd(Z8[1:length(monomials(x, [0, 1, 2, 3, 4])), :])
# r > 7
~, S = svd(Z8)

# Check for 2-shift-invariance:
~, S = svd(Z9[1:length(monomials(x, [0, 1, 2, 3])), :])
# r = 7 --> We have 2-shift-invariance! 

# ---------- Moment matrix approach (real roots only) ------------

# import SCS
# solver = SCS.Optimizer
include("solvers.jl")
big_clarabel = clarabel_optimizer(T = BigFloat)

# Start from Z_7:
d_mac = 7
Z = nullspace(macaulay(sys, d_mac, sparse_columns = false))
d_m = 3 # d_m = d_mom / 2
ν = moment_matrix(Z, big_clarabel, d_m, T = BigFloat)
# Clarabel: nearly optimal solution (so solver had some issues...)

# Manually inspect rank:
M = value_matrix(ν)
M, S = svd(M)
round.(log10.(S))

# Given the rank r, check for solutions:
r = 6
res = atomic_measure(ν, FixedRank(r), ShiftNullspace())
# The real-valued solution is identified! 

# ---------- Let's see whether a lower degree Macaulay suffices ----- 
d_mac = 6
Z = nullspace(macaulay(sys, d_mac, sparse_columns = false))
d_m = 3 # d_m = d_mom / 2
ν = moment_matrix(Z, big_clarabel, d_m, T = BigFloat)

# Manually inspect rank:
M = value_matrix(ν)
M, S = svd(M)
round.(log10.(S))

# Given the rank r, check for solutions:
r = 1
res = atomic_measure(ν, FixedRank(r), ShiftNullspace())
# The real-valued solution is identified! 

# d_mac = 6, d_mom = 3, d_tr = 3, r = 1: correct solution, but rank decision is non-trivial.

# Let's inspect flatness via truncations:
ν_tr = truncate(ν, 1)
M, S = svd(value_matrix(ν_tr))
round.(log10.(S))
r = 1
res = atomic_measure(ν_tr, FixedRank(r), ShiftNullspace())

# d_mac = 6, d_m = 3, d_trunc = 0, r = 1: no solution.
# d_mac = 6, d_m = 3, d_trunc = 1, r = 1: correct solution, easy rank decision.
# d_mac = 6, d_m = 3, d_trunc = 2, r = 1: correct solution, easy rank decision.

# Remember: in manuscript: d_mom = d_m * 2, d_tr = d_trunc * 2

# ----------- Old depreciated code ---------------

# This solves the feasibility problem - find moments that satisfy the given system
M = moment_matrix(sys, big_clarabel, 3)
# You provide t which determines the Macaulay matrix used to compute K_t (equivalent to null(Mac_t))
# Then, the function automatically computes M_d where d = roundDown(t/2)
# 3 / 2 = 1.5 rounded down = 1, so only degree 1 monomials in this moment matrix

# Try to find an atomic measure that corresponds to these moments
# Performs rank checks on principal submatrices, if flat --> setup eigenvalue problem to find atoms.
# Atoms of the measure corresponds to points v in V_c(I).
# ShiftNullspace() option computes atoms via shift-invariance of Macaulay null-space
atomic_measure(M, 1e-4, ShiftNullspace())
