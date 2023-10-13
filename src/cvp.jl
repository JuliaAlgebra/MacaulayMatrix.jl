import ManifoldsBase
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

function ManifoldsBase.retract!(
    M::Ellipsoid,
    q,
    p,
    x,
    t::Number,
    _::Manopt.AbstractRetractionMethod,
)
    return ManifoldsBase.project!(M, q, p + t * x)
end

# Not sure if that is correct
function ManifoldsBase.norm(
    ::Ellipsoid,
    p,
    x,
)
    return sqrt(x'x)
end

# Projection of `y` to the tangent space at `x`
function ManifoldsBase.project!(
    M::Ellipsoid,
    ret,
    x,
    y,
)
    z = M.A' * (M.A * x)
    ret .= y - z * (z'y) / (z'z)
    @assert isapprox(ret' * z, 0, atol = 1e-8)
    return ret
end

# Projection of `x` to the manifold
function ManifoldsBase.project!(
    M::Ellipsoid,
    ret,
    x,
)
    ret .= x ./ sqrt(x' * M.A' * M.A * x)
    @assert ret' * M.A' * M.A * ret ≈ 1
    return ret
end

function ManifoldsBase.zero_vector!(
    ::Ellipsoid,
    z,
    ::AbstractVector,
)
    z .= 0
    return z
end

function eliminate_indices(
    A::AbstractMatrix,
    I::AbstractVector{<:Integer};
    descent_state_type = Manopt.GradientDescentState,
    start = rand(size(A, 1)),
    stepsize = Manopt.ConstantStepsize(0.1),
    kws...,
)
    J = setdiff(axes(A, 2), I)
    AI = A[:, I]
    AJ = A[:, J]
    function eval_f_cb(M, x)
        val = norm(AI' * x)^2 / 2
        return val
    end
    function eval_grad_f_cb(M, x)
        grad = ManifoldsBase.project(M, x, AI * (AI' * x))
        return grad
    end
    manifold = Ellipsoid(AJ')
    mgo = Manopt.ManifoldGradientObjective(eval_f_cb, eval_grad_f_cb)
    dmgo = Manopt.decorate_objective!(manifold, mgo)
    problem = Manopt.DefaultManoptProblem(manifold, dmgo)
    s = descent_state_type(manifold, start; stepsize, kws...)
    state = Manopt.decorate_state!(s)
    Manopt.solve!(problem, state)
    obj = Manopt.get_objective(problem)
    return A' * Manopt.get_solver_return(obj, state)
end

function cvp(p, maxdegree; kws...)
    σ = MP.similar_variable(p, :σ)
    vars = MP.variables(p)
    grad = MP.differentiate(p, vars)
    eq = [σ - p; grad]
    M = macaulay(eq, maxdegree)
    N = eachindex(M.column_basis.monomials)
    non_pure_σ = findall(N) do i
        mono = M.column_basis.monomials[i]
        MP.degree(mono) != MP.degree(mono, σ)
    end
    J = setdiff(N, non_pure_σ)
    v = eliminate_indices(Matrix(SparseArrays.sparse(M)), non_pure_σ)
    return MP.polynomial(v[J], M.column_basis.monomials[J]), LinearAlgebra.norm(v[non_pure_σ])
end
