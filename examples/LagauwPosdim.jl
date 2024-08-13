using LinearAlgebra
using TypedPolynomials
using MacaulayMatrix
using JuMP
using MultivariateMoments
using DynamicPolynomials

# ----------------------------- Example 1 --------------------------------
# Small scale example with positive dimensional solution set over C. 

TypedPolynomials.@polyvar x y
sys_posdim = [
    x^2 + y^2
]

# Macaulay matrix method will never retrieve V_r(I), since V_c(I) is posdim:
sols = solve_system(sys_posdim, column_maxdegree=15, print_level=3)

# Nullity of Macaulay matrix increases with constant steps of 2
# this shows that the affine solution set is 1-dimensional.

# Let's try the moment-matrix approach:
d_max = 4
Z = nullspace(macaulay(sys_posdim, d_max))

# Or put them in a loop:
import SCS
solver = SCS.Optimizer
for i = 1:floor(d_max / 2)
    println("\n \n \n --------- M_$i --------- \n")
    M = moment_matrix(Z, solver, i)
    res = atomic_measure(M, 1e-4, ShiftNullspace())
    if length(res.atoms) >= 1
        println("-- Retrieved solutions: ")
        for j in eachindex(res.atoms)
            println(res.atoms[j].center)
        end
    end
end

# We are able to retrieve the real-valued solution.
# As v_r(I) = (0,0), the moment matrix is almost zero.


# ----------------------------- Example 2 --------------------------------