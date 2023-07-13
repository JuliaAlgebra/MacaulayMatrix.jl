module Macaulay

import LinearAlgebra, SparseArrays
import MultivariatePolynomials as MP
import MultivariateBases as MB
import MathOptInterface as MOI
import JuMP
import MultivariateMoments as MM
import Printf

import CommonSolve: solve, solve!, init, step!

abstract type AbstractIteration end

struct ColumnDegreeIteration <: AbstractIteration end

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

function step!(s::Iterator, ::ColumnDegreeIteration)
    if s.solver.max_iter > 0 && size(s.stats, 1) >= s.solver.max_iter
        s.status = MOI.ITERATION_LIMIT
        return
    end
    mindegree = maximum(MP.maxdegree, s.matrix.polynomials)
    added = 0
    for deg in mindegree:s.solver.column_maxdegree
        if deg == mindegree
            min = minimum(MP.maxdegree, s.matrix.polynomials)
            added = fill_column_maxdegrees!(s.matrix, min:deg)
        else
            added = fill_column_maxdegree!(s.matrix, deg)
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

include("hankel.jl")
include("monomial_generators.jl")

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
