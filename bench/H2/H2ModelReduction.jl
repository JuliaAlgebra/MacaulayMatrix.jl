# H2 globally optimal model reduction using Walsh's theorem
# Theory available on: https://ftp.esat.kuleuven.be/pub/stadius/slagauw/23-47.pdf
#
# (c) Sibren Lagauw, 2023.

using DynamicPolynomials
using HomotopyContinuation

include("auxiliaryFunctions.jl")

# Select higher-order system and degree of approximation
tf = selectProblem("F3")
m = 3

# Create system of equations:
system, vars = composeSystemWalsh(tf, m)

# Solve using Homotopy:
P = System(system, variables=vars)
@timev res = solve(P, only_non_zero=true, show_progress=true, threading=true, compile=true, catch_interrupt=true) # @timev (verbose option)
x = real_solutions(res)
T, validSols = analyseSolutions(x, tf, m, vars)

# Inspect valid stationary points (real-valued and stable):
T[validSols]