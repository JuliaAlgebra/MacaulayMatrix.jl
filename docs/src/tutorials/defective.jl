# # Defective univariate example
#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/defective.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/defective.ipynb)

using Test #src
using LinearAlgebra
using TypedPolynomials
using MacaulayMatrix
using MultivariateMoments
import Clarabel

# Consider the following system with a root of multiplicity 4:

@polyvar x
system = [(x - 4)^4]

# The Macaulay framework finds the root with degree 4:

solve_system(system, column_maxdegree = 4)

# Let's find the max rank PSD hankel from from the degree 4
# Macaulay nullspace:

ν = moment_matrix(system, Clarabel.Optimizer, 4)
nothing #hide

# We find the following PSD hankel

ν

# Looking at its singular values, the rank seems to be 1:

svd(value_matrix(ν)).S

# The rank of the truncation is also 1:

svd(value_matrix(truncate(ν, 1))).S

# We conclude that the PSD hankel is flat.
# The root can be obtained using:

η = #src
atomic_measure(ν4, FixedRank(1))

@test length(η.atoms) == 1 #src
@test η.atoms[1].center ≈ [4] rtol=1e-3 #src

# The solution is also apparent from the equation in the nullspace of the moment matrix:

nullspace(ν4, ShiftNullspace(), FixedRank(1)).polynomials
