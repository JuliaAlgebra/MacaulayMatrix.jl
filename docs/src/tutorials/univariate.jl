# # Univariate

using LinearAlgebra
using TypedPolynomials
using Macaulay
using JuMP
using MultivariateMoments

# Consider the following example:

@polyvar x
p = 3x^4 + 8x^3 - 6x^2 + 24x + 1
q = differentiate(p, x)

# Let's solve that system:

sols = solve_system([q], 4)
nothing # hide

# The real roots are:

sols

# With PSD Hankel matrix trick no complex root!

import CSDP
solver = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)
psd_hankel([q], solver, 4)

# What happened there ?
# First, we computed the Macaulay matrix

M, monos = macaulay_monomials([q], 4)
M

# We then get the nullspace

Z = nullspace(Matrix(M))

# The PSD hankel matrix we find is:

M = moment_matrix(Z, solver, 2, monos)

# From which we get:

extractatoms(M, 1e-6)
