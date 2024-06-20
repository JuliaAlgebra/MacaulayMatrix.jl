# H2 globally optimal model reduction using Walsh's theorem
# Theory available on: https://ftp.esat.kuleuven.be/pub/stadius/slagauw/23-47.pdf
#
# (c) Sibren Lagauw, 2023.

using MacaulayMatrix

# Create system of equations:
m = 3
system, vars = H2.composeSystemWalsh(H2.F3(), m)

# Solve using Homotopy:
using HomotopyContinuation
P = System(system, variables=vars)
@timev res = solve(P, only_non_zero=true, show_progress=true, threading=true, compile=true, catch_interrupt=true) # @timev (verbose option)
x = real_solutions(res)
T, validSols = H2.analyseSolutions(x, tf, m, vars)

# Inspect valid stationary points (real-valued and stable):
T[validSols]
