# # Univariate

using LinearAlgebra
using DynamicPolynomials
using MultivariateMoments
using JuMP
using MacaulayMatrix

# Consider the following example:

@polyvar x
p = 3x^4 + 8x^3 - 6x^2 + 24x + 1
q = differentiate(p, x)

# Let's solve that system:

sols = solve_system([q], column_maxdegree = 4)
nothing # hide

# The real roots are:

sols

# With PSD Hankel matrix trick no complex root!

import SCS
solver = SCS.Optimizer
psd_hankel([q], solver, 4)

# What happened there ?
# First, we computed the MacaulayMatrix matrix nullspace


Z = nullspace(macaulay([q], 4))

# The PSD hankel matrix we find is:

M = moment_matrix(Z, solver, 2)

# From which we get:

atomic_measure(M, 1e-6)
