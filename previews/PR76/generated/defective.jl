using LinearAlgebra
using TypedPolynomials
using MacaulayMatrix
using MultivariateMoments
import Clarabel

@polyvar x
system = [(x - 4)^4]

solve_system(system, column_maxdegree = 4)

ν = moment_matrix(system, Clarabel.Optimizer, 4)
nothing #hide

ν

svd(value_matrix(ν)).S

svd(value_matrix(truncate(ν, 1))).S

atomic_measure(ν, FixedRank(1))

nullspace(ν, ShiftNullspace(), FixedRank(1)).polynomials

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
