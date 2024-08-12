using LinearAlgebra
using TypedPolynomials
using MacaulayMatrix
using JuMP
using MultivariateMoments

@polyvar x[1:4]
system = [
    x[1] + x[2] - 1,
    x[1] * x[3] + x[2] * x[4],
    x[1] * x[3]^2 + x[2] * x[4]^2 - 2/3,
    x[1] * x[3]^3 + x[2] * x[4]^3,
]

solve_system(system, column_maxdegree = 8)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
