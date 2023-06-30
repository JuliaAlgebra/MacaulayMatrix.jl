module Macaulay

import LinearAlgebra, SparseArrays
import MultivariatePolynomials as MP
import MultivariateBases as MB
import JuMP
import MultivariateMoments as MM
import Printf

function realization_observability(U::AbstractMatrix)
    if true
        #shift = size(U, 2)
        shift = 1
        U_ = U[1:(end-shift), :]
        _U = U[(shift + 1):end, :]
        A = LinearAlgebra.pinv(U_) * _U
    else
        r = size(U, 2)
        U_ = U[1:r, :]
        _U = U[2:(r+1), :]
        A = LinearAlgebra.pinv(U_) * _U
    end
    return [[x] for x in LinearAlgebra.eigen(A).values]
end

function old_realization_hankel(H::LinearAlgebra.Symmetric)
    ldlt = MM.low_rank_ldlt(
        H,
        MM.SVDLDLT(),
        MM.LeadingRelativeRankTol(1e-6),
    )
    return realization_observability(ldlt.L)
end

function realization_hankel(M::MM.MomentMatrix)
    η = MM.atomic_measure(M, 1e-4)
    if isnothing(η)
        return
    else
        return [η.atoms[i].center for i in eachindex(η.atoms)]
    end
end

function realization_hankel(H::LinearAlgebra.Symmetric, monos)
    return realization_hankel(MM.MomentMatrix(H, monos))
end

function MM.moment_matrix(Z::MM.MacaulayNullspace, solver, d; print_level=1)
    model = JuMP.Model(solver)
    r = size(Z.matrix, 2)
    JuMP.@variable(model, b[1:r])
    JuMP.@constraint(model, sum(b) == 1)
    Zb = Z.matrix * b
    monos = Z.basis.monomials
    gram_monos = MP.monomials(MP.variables(monos), 0:d)
    H = LinearAlgebra.Symmetric([
        Zb[findfirst(isequal(gram_monos[i] * gram_monos[j]), monos)]
        for i in eachindex(gram_monos), j in eachindex(gram_monos)
    ])
    JuMP.@constraint(model, H in JuMP.PSDCone())
    JuMP.optimize!(model)
    if print_level >= 1
        @info("Terminated with $(JuMP.termination_status(model)) ($(JuMP.raw_status(model))) in $(JuMP.solve_time(model)) seconds.")
    end
    if JuMP.termination_status(model) == JuMP.MOI.INFEASIBLE
        return
    elseif JuMP.termination_status(model) == JuMP.MOI.ALMOST_OPTIMAL
        @warn(string(JuMP.solution_summary(model)))
    elseif JuMP.termination_status(model) != JuMP.MOI.OPTIMAL
        error(string(JuMP.solution_summary(model)))
    end
    H = LinearAlgebra.Symmetric(JuMP.value.(H))
    return MM.MomentMatrix(H, gram_monos)
end

function MM.moment_matrix(polynomials::AbstractVector{<:MP.AbstractPolynomialLike{T}}, solver, maxdegree) where {T}
    Z = macaulay_nullspace(polynomials, maxdegree)
    return MM.moment_matrix(Z, solver, div(maxdegree, 2))
end

function psd_hankel(args...)
    H = MM.moment_matrix(args...)
    if H === nothing
        return nothing
    end
    return realization_hankel(H)
end

macaulay(polynomials, maxdegree) = first(macaulay_monomials(polynomials, maxdegree))

function macaulay_monomials(polynomials::AbstractVector{<:MP.AbstractPolynomialLike{T}}, maxdegree) where {T}
    vars = MP.variables(polynomials)
    monos = MP.monomials(vars, 0:maxdegree)
    column = Dict(monos[i] => i for i in eachindex(monos))
    row = 0
    I = Int[]
    J = Int[]
    K = T[]
    for leading_mono in monos
        for p in polynomials
            lm = MP.leading_monomial(p)
            if MP.divides(lm, leading_mono)
                factor = MP.div_multiple(leading_mono, lm)
                row += 1
                for t in MP.terms(p)
                    push!(I, row)
                    push!(J, column[factor * MP.monomial(t)])
                    push!(K, MP.coefficient(t))
                end
            end
        end
    end
    return SparseArrays.sparse(I, J, K, row, length(monos)), monos
end

function _nullspace(
    M::Matrix,
    # This corresponds to the default of `LinearAlgebra.nullspace`
    rank_check=MM.LeadingRelativeRankTol(min(size(M)...) * eps(real(float(oneunit(eltype(M)))))),
)
    m, n = size(M)
    if iszero(m) || iszero(n)
        Z = Matrix{LinearAlgebra.eigtype(eltype(A))}(LinearAlgebra.I, n, n)
        accuracy = zero(T)
    else
        SVD = LinearAlgebra.svd(M; full=true)
        r = MM.rank_from_singular_values(SVD.S, rank_check)
        Z = (SVD.Vt[(r+1):end,:])'
        accuracy = MM.accuracy(SVD.S, r, rank_check)
    end
    return Z, accuracy
end

_nullspace(M, args...) = _nullspace(Matrix(M), args...)

function macaulay_nullspace(polynomials::AbstractVector{<:MP.AbstractPolynomialLike{T}}, maxdegree, args...) where {T}
    Δt = @elapsed begin
        M, monos = macaulay_monomials(polynomials, maxdegree)
        Z, accuracy = _nullspace(M, args...)
    end
    @info("Nullspace of degree $maxdegree of dimensions $(size(Z)) computed in $Δt seconds.")
    return MM.MacaulayNullspace(Z, MB.MonomialBasis(monos), accuracy)
end

# Inspired from `macaulaylab.net/Code/solvesystemnullspace.m`
function solve_system(polynomials::AbstractVector{<:MP.AbstractPolynomialLike{T}}, maxdegree, args...; print_level=1) where {T}
    mindegree = maximum(MP.maxdegree, polynomials)
    nullities = zeros(Int, maxdegree)
    Z = nothing
    Printf.@printf("\t | degree \t | nullity \t | increase \t |\n")
    Printf.@printf("\t |-----------------------------------------------|\n")
    for d in mindegree:maxdegree
        Z = macaulay_nullspace(polynomials, d, args...)
        nullities[d] = size(Z.matrix, 2)
        change = nullities[d] - nullities[d - 1]
        if print_level >= 1
            Printf.@printf(
                "\t | %d \t \t | %d \t \t | %d \t \t |\n",
                d,
                nullities[d],
                change,
            )
        end
        if d > mindegree && iszero(change)
            sols = MM.solve(Z, MM.ShiftNullspace())
            if !isnothing(sols)
                return sols.elements
            end
        end
    end
    return
end

import SemialgebraicSets as SS
struct Solver <: SS.AbstractAlgebraicSolver
    maxdegree::Int
    print_level::Int
end
Solver(maxdegree) = Solver(maxdegree, 1)

SS.default_gröbner_basis_algorithm(::Any, ::Solver) = SS.NoAlgorithm()

SS.promote_for(::Type{T}, ::Type{Solver}) where {T} = float(T)

function SS.solve(V::SS.AbstractAlgebraicSet, solver::Solver)
    return solve_system(V.polynomials, solver.maxdegree; print_level = solver.print_level)
end


# Taken from JuMP.jl
# This package exports everything except internal symbols, which are defined as those
# whose name starts with an underscore. Macros whose names start with
# underscores are internal as well. If you don't want all of these symbols
# in your environment, then use `import JuMP` instead of `using JuMP`.

const _EXCLUDE_SYMBOLS = [Symbol(@__MODULE__), :eval, :include]

for sym in names(@__MODULE__; all = true)
    sym_string = string(sym)
    if sym in _EXCLUDE_SYMBOLS ||
       startswith(sym_string, "_") ||
       startswith(sym_string, "@_")
        continue
    end
    if !(
        Base.isidentifier(sym) ||
        (startswith(sym_string, "@") && Base.isidentifier(sym_string[2:end]))
    )
        continue
    end
    @eval export $sym
end

include("H2/H2.jl")
export H2

end # module Macaulay
