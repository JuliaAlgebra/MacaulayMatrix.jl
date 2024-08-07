# From Example 1 of "Exact characterization of the global optima of least squares realization of autonomous LTI models as a multiparameter eigenvalue problem"

using TypedPolynomials
@polyvar a[1:1] g[1:2]
system = [
    4a[1] - 2a[1] * g[1] - a[1]^2 * g[2] - a[1]^3 * g[1] + 3,
    3a[1] - g[1] - 2a[1] * g[2] - 2a[1]^2 * g[1] - a[1]^3 * g[2] + 2,
    2a[1] - g[2] - a[1] * g[1] - 2a[1]^2 * g[2] + 1,
]

uni = [2a[1]^7 - 5a[1]^6 + 12a[1]^5 - a[1]^4 + 6a[1]^3 + 3a[1]^2 + 3]

using MacaulayMatrix
using MultivariateMoments
solve_system(uni, column_maxdegree = 7)

solve_system(system, column_maxdegree = 7)
it = Iterator(system, Solver())
step!(it)
step!(it)
step!(it)
step!(it)
using LinearAlgebra
null = nullspace(it.matrix)
null.matrix
sols = solve(null, ShiftNullspace())

using Clarabel
ν6 = moment_matrix(system, Clarabel.Optimizer, 6)
for i in 1:20
    @show atomic_measure(ν6, FixedRank(i))
end

using Clarabel
ν7 = moment_matrix(system, Clarabel.Optimizer, 7)
atomic_measure(ν7, FixedRank(5))
real_null = MacaulayNullspace(ν, FixedRank(5))
sols = solve(real_null, ShiftNullspace())
border, solver = border_basis_and_solver(real_null, ShiftNullspace())
@edit solve(border, solver)

import SCS
ν8 = moment_matrix(system, SCS.Optimizer, 8)
for r in 1:size(ν8.Q, 1)
    @show r
    @show atomic_measure(ν8, FixedRank(r))
end

using MosekTools
ν8 = moment_matrix(system, Mosek.Optimizer, 8)
@show atomic_measure(ν8, LeadingRelativeRankTol(1e-1))
real_null8 = MacaulayNullspace(ν8, FixedRank(5))
sols = solve(real_null8, ShiftNullspace())
function truncate(ν, d)
    μ = MultivariateMoments.measure(ν)
    monos = monomials(variables(μ), 0:d)
    MultivariateMoments.hankel(μ, monos, monos)
end
using MultivariateBases
function mtruncate(ν, d)
    Q = truncate(ν, d)
    μ = MultivariateMoments.measure(ν)
    monos = monomials(variables(μ), 0:d)
    return MomentMatrix(Q, monos)
end
truncated_rank(ν, d) = rank(truncate(ν, d))
border8, solver = border_basis_and_solver(real_null8, ShiftNullspace())
rank(truncate(ν8, 1), rtol=1e-1)
rank(truncate(ν8, 2), rtol=1e-1)
rank(truncate(ν8, 3), rtol=1e-1)
rank(truncate(ν8, 4), rtol=1e-1)
rank(Matrix(ν8.Q), rtol=1e-1)
ν8_3 = mtruncate(ν8, 3)
atomic_measure(ν8_3, FixedRank(1))
atomic_measure(ν8_3, FixedRank(2))
atomic_measure(ν8_3, FixedRank(3))
atomic_measure(ν8_3, LeadingRelativeRankTol(1e-3))
svd(Matrix(ν8_3.Q))
ν8_3
svd(truncate(ν8, 3))

ν10 = moment_matrix(system, Mosek.Optimizer, 10)
rank(truncate(ν10, 1), rtol=1e-1)
rank(truncate(ν10, 2), rtol=1e-1)
rank(truncate(ν10, 3), rtol=1e-1)
rank(truncate(ν10, 4), rtol=1e-1)
rank(truncate(ν10, 5), rtol=1e-1)
ν12 = moment_matrix(system, Mosek.Optimizer, 12)
rank(truncate(ν12, 1), rtol=1e-1)
rank(truncate(ν12, 2), rtol=1e-1)
rank(truncate(ν12, 3), rtol=1e-1)
rank(truncate(ν12, 4), rtol=1e-1)
rank(truncate(ν12, 5), rtol=1e-1)
rank(truncate(ν12, 6), rtol=1e-1)

