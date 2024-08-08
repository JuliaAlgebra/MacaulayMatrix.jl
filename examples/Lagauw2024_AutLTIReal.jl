using LinearAlgebra
using TypedPolynomials
using MacaulayMatrix
using JuMP
using MultivariateMoments
using HomotopyContinuation


# Walsh MEP from (Lagauw, Vanpoucke, De Moor, 2024)
function composeSystemWalsh(y::Vector{Float64}, n::Int)
    @var a[1:n]
    aa = [1; a]

    N = length(y)
    T = zeros(Expression, N - n, N)
    for i in 1:N-n
        T[i, i:i-1+length(aa)] = reverse(aa)
    end
    S = zeros(Expression, N - 2n, N - n)
    for i in 1:N-2n
        S[i, i:i-1+length(aa)] = reverse(aa)
    end

    @var g[1:N-2n]

    system = Vector{Expression}(undef, N - n)
    system[1:N-n] = T * y - T * transpose(T) * transpose(S) * g

    return system
end


n = 1;
y = [4.0, 3, 2, 1];
# D = 3 (degree of polynomials is at most cubic)
sys = composeSystemWalsh(y, n)

#@var a g[1:2]
#system = [
#    3 + 4 * a - (a * (g[1] + g[2] * a) + g[1] * a + g[1] * a^3)
#    2 + 3 * a - (g[1] + a * (a * (g[1] + g[2] * a) + g[1] * a) + 2 * g[2] * a)
#    1 + 2 * a - (g[2] + a * (g[1] + 2 * g[2] * a))
#]


# ---------- Standard Macaulay Null space ------------

# Degree 7 Macaulay matrix suffices to retrieve V_r(I), by selection from V_c(I):
sols = solve_system(system, column_maxdegree=7, print_level=3)

# Step by step analysis:
solver = Iterator(system, MacaulayMatrix.Solver())

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

# using Plots
# plot(saturated_dependence(solver))

# ---------- Moment matrix approach (real roots only) ------------

import SCS
solver = SCS.Optimizer


# This solves the feasibility problem - find moments that satisfy the given system
M = moment_matrix(system, solver, 3)
# You provide t which determines the Macaulay matrix used to compute K_t (equivalent to null(Mac_t))
# Then, the function automatically computes M_d where d = roundDown(t/2)
# 3 / 2 = 1.5 rounded down = 1, so only degree 1 monomials in this moment matrix

# Try to find an atomic measure that corresponds to these moments
# Performs rank checks on principal submatrices, if flat --> setup eigenvalue problem to find atoms.
# Atoms of the measure corresponds to points v in V_c(I).
# ShiftNullspace() option computes atoms via shift-invariance of Macaulay null-space
atomic_measure(M, 1e-4, ShiftNullspace())

# t = D = 3 did not suffice...

# t = 6?
M = moment_matrix(system, solver, 6)
atomic_measure(M, 1e-4, ShiftNullspace())

# You need K_t with t = 7, however, M_3(y) with y \in K_7 works!
# This function automatically computes M_d where d = roundDown(t/2)
M = moment_matrix(system, solver, 7)
atomic_measure(M, 1e-4, ShiftNullspace())

# You can also first compute the nullspace, than provide it to compute a moment_matrix
Z = nullspace(macaulay(system, 7))
M = moment_matrix(Z, solver, 3)
atomic_measure(M, 1e-4, ShiftNullspace())


# TODO: play with rank options: FixedRank(), ... See most recent push Benoit. 

# Idea: Use rewriting family of real radical ideal to see if we can find equations that eliminate complex solutions and add them to the system of polynomial equations! This equation might generalize to other instances of y!

