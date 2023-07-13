module Macaulay

import LinearAlgebra, SparseArrays
import MultivariatePolynomials as MP
import MultivariateBases as MB
import MathOptInterface as MOI
import JuMP
import MultivariateMoments as MM
import Printf

import CommonSolve: solve, solve!, init, step!

"""
    add_monomial_ideal_generators!(generators, standard_monomials, fixed, vars)

Given the standard monomials of a zero-dimensional monomial ideal, returns the
generators of this monomial ideal of the form `fixed * m` where `m` is a
monomial in `vars` and `vars * m` is not a multiple of any of the monomials in
`generators`.
"""
function add_monomial_ideal_generators!(generators, standard_monomials::MB.MonomialBasis, fixed, vars)
    if any(g -> MP.divides(g, fixed), generators)
        return
    end
    if isempty(vars)
        push!(generators, fixed)
    else
        var = last(vars)
        maxdeg = maximum(
            (
                MP.degree(mono, var)
                for mono in standard_monomials.monomials
                if all(v -> v in vars || MP.degree(mono, v) == MP.degree(fixed, v), MP.variables(fixed))
            );
            init = -1,
        )
        maxdeg += 1
        if length(vars) == 1
            degs = maxdeg:maxdeg
        else
            degs = 0:maxdeg
        end
        other_vars = vars[1:end-1]
        for d in degs
            add_monomial_ideal_generators!(generators, standard_monomials, fixed * var^d, other_vars)
        end
    end
end

"""
    monomial_ideal_generators(standard_monomials)

Given the standard monomials of a zero-dimensional monomial ideal, returns the
generators of this monomial ideal.
"""
function monomial_ideal_generators(standard_monomials::MB.MonomialBasis)
    generators = similar(standard_monomials.monomials, 0)
    fixed = MP.constant_monomial(eltype(standard_monomials.monomials))
    add_monomial_ideal_generators!(generators, standard_monomials, fixed, MP.variables(standard_monomials))
    return MB.MonomialBasis(generators)
end

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

function MM.moment_matrix(null::MM.MacaulayNullspace, solver, d; print_level=1)
    # TODO Newton polytope
    vars = MP.variables(null.basis.monomials)
    monos = MP.monomials(vars, 0:2d)
    model = JuMP.Model(solver)
    r = size(null.matrix, 2)
    # Number of roots at infinity that were left out by not adding these as columns
    null = null[[mono for mono in monos if mono in null.basis.monomials]]
    num_inf = length(monos) - size(null.matrix, 1)
    JuMP.@variable(model, b[1:(r + num_inf)])
    JuMP.@constraint(model, sum(b) == 1)
    inf_idx = 0
    Zb = map(monos) do mono
        idx = MM._monomial_index(null.basis.monomials, mono)
        if isnothing(idx)
            inf_idx += 1
            return convert(JuMP.AffExpr, b[r + inf_idx])
        else
            return null.matrix[idx, :]' * b[1:r]
        end
    end
    @assert inf_idx == num_inf
    gram_monos = MP.monomials(vars, 0:d)
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
    Z = LinearAlgebra.nullspace(macaulay(polynomials, maxdegree))
    return MM.moment_matrix(Z, solver, div(maxdegree, 2))
end

function psd_hankel(args...)
    H = MM.moment_matrix(args...)
    if H === nothing
        return nothing
    end
    return realization_hankel(H)
end

abstract type AbstractIteration end

Base.@kwdef struct ColumnDegreeIteration <: AbstractIteration
    sparse_columns::Bool = true
end

import SemialgebraicSets as SS
# `@kwdef` is not exported in Julia v1.6
Base.@kwdef mutable struct Solver <: SS.AbstractAlgebraicSolver
    default_iteration::AbstractIteration = ColumnDegreeIteration()
    column_maxdegree::Int = 0
    print_level::Int = 1
    max_iter::Int = 10
    rank_check::Union{Nothing,MM.RankCheck} = nothing
end

SS.default_gröbner_basis_algorithm(::Any, ::Solver) = SS.NoAlgorithm()

SS.promote_for(::Type{T}, ::Type{Solver}) where {T} = float(T)

import DataFrames

include("matrix.jl")

mutable struct Iterator{
    T,
    P<:MP.AbstractPolynomialLike,
    V<:AbstractVector{P},
    B,
    U,
}
    matrix::MacaulayMatrix{T,P,V,B}
    standard_monomials::Union{Nothing,B}
    border_monomials::Union{Nothing,B}
    solutions::Union{Nothing,Vector{Vector{U}}}
    status::MOI.TerminationStatusCode
    stats::DataFrames.DataFrame
    solver::Solver
    function Iterator(
        matrix::MacaulayMatrix{T,P,V,B},
        solver::Solver,
    ) where {T,P,V,B}
        U = SS.promote_for(T, Solver)
        return new{T,P,V,B,U}(
            matrix,
            nothing,
            nothing,
            nothing,
            MOI.OPTIMIZE_NOT_CALLED,
            DataFrames.DataFrame(nullity = Int[], num_rows = Int[], num_cols = Int[]),
            solver,
        )
    end
end

function Iterator(
    polynomials::AbstractVector{<:MP.AbstractPolynomialLike},
    solver::Solver,
)
    return Iterator(MacaulayMatrix(polynomials), solver)
end

function Base.show(io::IO, s::Iterator)
    println(io, "Macaulay matrix solver. Last iteration considered:")
    show(io, s.matrix)
    println(io, "Current status is $(s.status)")
    if !isnothing(s.solutions)
        println(io, "Found $(length(s.solutions)) solutions:")
        for sol in s.solutions
            println(io, "  ", sol)
        end
    end
    println(io, "History of iterations:")
    show(io, s.stats)
end

init(V::SS.AbstractAlgebraicSet, solver::Solver) = Iterator(SS.equalities(V), solver)

function solve!(s::Iterator)
    while s.status == MOI.OPTIMIZE_NOT_CALLED
        step!(s)
    end
    return s.solutions
end

# Inspired from `macaulaylab.net/Code/solvesystemnullspace.m`
function solve_system(polynomials::AbstractVector{<:MP.AbstractPolynomialLike}; kws...)
    return solve(SS.algebraic_set(polynomials, Solver(; kws...)))
end

step!(s::Iterator) = step!(s, s.solver.default_iteration)

function step!(s::Iterator, it::ColumnDegreeIteration)
    if s.solver.max_iter > 0 && size(s.stats, 1) >= s.solver.max_iter
        s.status = MOI.ITERATION_LIMIT
        return
    end
    mindegree = maximum(MP.maxdegree, s.matrix.polynomials)
    added = 0
    for deg in mindegree:s.solver.column_maxdegree
        if deg == mindegree
            min = minimum(MP.maxdegree, s.matrix.polynomials)
            added = fill_column_maxdegrees!(s.matrix, min:deg; sparse_columns = it.sparse_columns)
        else
            added = fill_column_maxdegree!(s.matrix, deg; sparse_columns = it.sparse_columns)
        end
        if !iszero(added)
            break
        end
    end
    if iszero(added)
        s.status = MOI.OTHER_LIMIT
        return
    end
    if s.solver.print_level >= 3
        display(s)
    end
    args = isnothing(s.solver.rank_check) ? tuple() : (s.solver.rank_check,)
    Z = LinearAlgebra.nullspace(s.matrix, args...)
    nullity = size(Z.matrix, 2)
    if !isempty(s.stats.nullity) && nullity == s.stats.nullity[end]
        sols = MM.solve(Z, MM.ShiftNullspace())
        if !isnothing(sols)
            s.solutions = collect(sols)
            if isempty(s.solutions)
                s.status = MOI.INFEASIBLE
            else
                s.status = MOI.OPTIMAL
            end
        end
    end
    push!(s.stats, [nullity, size(s.matrix)...])
    return
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
