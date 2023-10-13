import ManifoldsBase
import Manifolds
import Manopt
export eliminate_indices

"""
    function Ellipsoid{T,M<:AbstractMatrix}
        M::M
    end

Represent the set `x'A'Ax = 1`.
"""
struct Ellipsoid{T,M<:AbstractMatrix{T}} <: ManifoldsBase.AbstractManifold{ManifoldsBase.ℝ}
    A::M
end

ManifoldsBase.manifold_dimension(e::Ellipsoid) = size(e.A, 2)

ManifoldsBase.inner(::Ellipsoid, p, x, y) = LinearAlgebra.dot(x, y)

# See Example 3.49, p. 40 of Boumal's book for the retraction on the Sphere
function ManifoldsBase.retract_project!(M::Ellipsoid, q, p, x, t::Number)
    return ManifoldsBase.project!(M, q, p + t * x)
end

# Projection of `y` to the tangent space at `x`
function ManifoldsBase.project!(M::Ellipsoid, ret, x, y)
    z = M.A' * (M.A * x)
    ret .= y - z * (z'y) / (z'z)
    @assert isapprox(ret' * z, 0, atol = 1e-8)
    return ret
end

# Projection of `x` to the manifold
function ManifoldsBase.project!(M::Ellipsoid, ret, x)
    ret .= x ./ sqrt(x' * M.A' * M.A * x)
    @assert ret' * M.A' * M.A * ret ≈ 1
    return ret
end

function ManifoldsBase.zero_vector!(::Ellipsoid, z, ::AbstractVector)
    z .= 0
    return z
end

import ForwardDiff

function _eliminate_indices(
    manifold::ManifoldsBase.AbstractManifold,
    A::AbstractMatrix,
    I::AbstractVector{<:Integer};
    partition,
    start = rand(size(A, 1)),
    solver = Manopt.GradientDescentState,
    kws...,
)
    AI = A[:, I]
    if solver === Manopt.LevenbergMarquardt
        function f(M, x)
            return AI' * x
        end
        function jacobian_f(M, p; basis_domain)
            #N = LinearAlgebra.nullspace(p.x[1]')
            #B = Matrix{Float64}(undef, size(AI, 1) - 1, size(AI, 2))
            #n = length(p.x[1])
            #B[1:(n-1), :] = N' * AI[1:n, :]
            #B[n:end, :] = AI[(n+1):end, :]
            #return B'
            X0 = zeros(ManifoldsBase.manifold_dimension(M))
            J = ForwardDiff.jacobian(
                x -> f(
                    M,
                    ManifoldsBase.exp(
                        M,
                        p,
                        ManifoldsBase.get_vector(M, p, x, basis_domain),
                    ),
                ),
                X0,
            )
            return J
        end
        x = Manopt.LevenbergMarquardt(manifold, f, jacobian_f, partition(start), length(I))
    else
        function eval_f_cb(M, x)
            val = LinearAlgebra.norm(AI' * x)^2 / 2
            return val
        end
        function eval_grad_f_cb(M, x)
            g = partition(AI * (AI' * x))
            grad = ManifoldsBase.project(M, x, g)
            return grad
        end
        mgo = Manopt.ManifoldGradientObjective(eval_f_cb, eval_grad_f_cb)
        dmgo = Manopt.decorate_objective!(manifold, mgo)
        problem = Manopt.DefaultManoptProblem(manifold, dmgo)
        s = solver(manifold, partition(start); kws...)
        state = Manopt.decorate_state!(s)
        Manopt.solve!(problem, state)
        obj = Manopt.get_objective(problem)
        x = Manopt.get_solver_return(obj, state)
    end
    return A' * x
end

function eliminate_indices(
    A::AbstractMatrix,
    I::AbstractVector{<:Integer};
    stepsize = Manopt.ConstantStepsize(0.1),
    retraction_method = ManifoldsBase.ProjectionRetraction(),
    kws...,
)
    J = setdiff(axes(A, 2), I)
    manifold = Ellipsoid(A[:, J]')
    return _eliminate_indices(
        manifold,
        A,
        I;
        stepsize,
        partition = identity,
        retraction_method,
        kws...,
    )
end

function eliminate_indices(
    A::AbstractMatrix,
    I::AbstractVector{<:Integer},
    num_nz_rows::Integer;
    kws...,
)
    num_z_rows = size(A, 1) - num_nz_rows
    partition(x) = Manifolds.ArrayPartition(x[1:num_nz_rows], x[(num_nz_rows+1):end])
    manifold = Manifolds.ProductManifold(
        Manifolds.Sphere(num_nz_rows - 1),
        Manifolds.Euclidean(num_z_rows),
    )
    return _eliminate_indices(manifold, A, I; partition, kws...)
end

function pure_power(mono, var)
    return MP.degree(mono) == MP.degree(mono, var)
end

function _zero_constant(p)
    return iszero(MP.coefficient(p, MP.constant_monomial(p)))
end

function pure_σ_rows(M::MacaulayMatrix, i::Int, σ::MP.AbstractVariable)
    σ_rows = Int[]
    other_rows = Int[]
    row = 0
    for (shift, statuses) in zip(M.row_shifts.monomials, M.shift_statuses)
        for (j, status) in enumerate(statuses)
            if status == INCLUDED || status == NOT_REDUNDANT
                row += 1
                if pure_power(shift, σ) && (j == i || !_zero_constant(M.polynomials[j]))
                    push!(σ_rows, row)
                else
                    push!(other_rows, row)
                end
            end
        end
    end
    return σ_rows, other_rows
end

function cvp(
    p,
    maxdegree;
    ellipsoid = false,
    rank_check = MM.LeadingRelativeRankTol(Base.rtoldefault(Float64)),
    kws...,
)
    σ = MP.similar_variable(p, :σ)
    vars = MP.variables(p)
    grad = MP.differentiate(p, vars)
    eq = [σ - p; grad]
    M = macaulay(eq, maxdegree)
    N = eachindex(M.column_basis.monomials)
    non_pure_σ = findall(N) do i
        mono = M.column_basis.monomials[i]
        return !pure_power(mono, σ)
    end
    J = setdiff(N, non_pure_σ)
    S = SparseArrays.sparse(M)
    if ellipsoid
        v = eliminate_indices(S, non_pure_σ; kws...)
    else
        σ_rows, other_rows = pure_σ_rows(M, 1, σ)
        reordered = [σ_rows; other_rows]
        S = S[reordered, :]
        n = length(σ_rows)
        # We now work in the manifold `Sphere(n-1) × R^(m-n)`
        if !isnothing(rank_check)
            A = S[1:n, J]
            F = LinearAlgebra.svd(Matrix(A), full = true)
            r = MM.rank_from_singular_values(F.S, rank_check)
            # We reformulate `x` as `U * y`.
            # Since `A' * U[:, (r+1):end] * y = 0`, we can reformulate as
            # `Sphere(r-1) × R^(m-r)`
            # TODO this destroys sparsity, we shouldn't actually do this product
            S = [F.U' * S[1:n, :]; S[(n+1):end, :]]
        else
            r = n
            F = nothing
        end
        v = eliminate_indices(S, non_pure_σ, r; kws...)
    end
    return MP.polynomial(v[J], M.column_basis.monomials[J]),
    LinearAlgebra.norm(v[non_pure_σ])
end
