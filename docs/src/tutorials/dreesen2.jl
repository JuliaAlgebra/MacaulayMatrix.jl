# # Getting rid of root at infinity

using LinearAlgebra
using TypedPolynomials
using Macaulay
using JuMP
import CSDP
using MultivariateMoments


# Consider the `dreesen2` example:

@polyvar x y z
system = [
    x^2 - x*y + z,
    2y^3 - 2x*y^2 - 3x * y,
    z^3 - x*y*z - 2,
]

# We first try to solve the system:

sols = solve_system(system, 6)
nothing # hide

# The real solutions are

sols

# With moment matrix of degree 3:

import CSDP
solver = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)
M = moment_matrix(system, solver, 3)
extractatoms(M, 1e-4, ShiftNullspace())

# With moment matrix of degree 4:

M = moment_matrix(system, solver, 4)
extractatoms(M, 1e-4, ShiftNullspace())

# With moment matrix of degree 5:

M = moment_matrix(system, solver, 5)
extractatoms(M, 1e-4, ShiftNullspace())
