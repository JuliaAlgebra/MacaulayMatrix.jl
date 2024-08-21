using LinearAlgebra
using DynamicPolynomials
using MultivariateMoments
using JuMP
using MacaulayMatrix

@polyvar x
p = 3x^4 + 8x^3 - 6x^2 + 24x + 1
q = differentiate(p, x)

sols = solve_system([q], column_maxdegree = 4)
nothing # hide

sols

import SCS
solver = SCS.Optimizer
psd_hankel([q], solver, 4)

Z = nullspace(macaulay([q], 4))

M = moment_matrix(Z, solver, 2)

atomic_measure(M, 1e-6)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
