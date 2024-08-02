using LinearAlgebra
using TypedPolynomials
using MacaulayMatrix
using JuMP
using MultivariateMoments

@polyvar x y z
system = [
    x^2 - x*y + z,
    2y^3 - 2x*y^2 - 3x * y,
    z^3 - x*y*z - 2,
]

sols = solve_system(system, column_maxdegree = 6)
nothing # hide

sols

solver = Iterator(system, MacaulayMatrix.Solver())
step!(solver)

solver

using Plots
plot(saturated_dependence(solver))

step!(solver)

plot(saturated_dependence(solver))

step!(solver)

plot(saturated_dependence(solver))

step!(solver)

plot(saturated_dependence(solver))

solver = Iterator(system, MacaulayMatrix.Solver())
step!(solver)

step!(solver, FirstStandardNonSaturated(10))

plot(saturated_dependence(solver))

step!(solver, FirstStandardNonSaturated(10))

plot(saturated_dependence(solver))

step!(solver, FirstStandardNonSaturated(10))

plot(saturated_dependence(solver))

step!(solver, FirstStandardNonSaturated(10))

plot(saturated_dependence(solver))

import SCS
solver = SCS.Optimizer
M = moment_matrix(system, solver, 3)
nothing # hide

atomic_measure(M, 1e-4, ShiftNullspace())

M = moment_matrix(system, solver, 4)
nothing # hide

atomic_measure(M, 1e-4, ShiftNullspace())

M = moment_matrix(system, solver, 5)
nothing # hide

atomic_measure(M, 1e-4, ShiftNullspace())

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
