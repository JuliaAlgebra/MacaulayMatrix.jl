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

import SemialgebraicSets as SS
# `@kwdef` is not exported in Julia v1.6
Base.@kwdef mutable struct Solver <: SS.AbstractAlgebraicSolver
    trim_to_border::Bool = false
    sparse_columns::Bool = true
    wait_for_gap::Bool = false
    default_iteration::Any = ColumnDegrees
    column_maxdegree::Int = 0
    print_level::Int = 1
    max_iter::Int = 10
    rank_check::Union{Nothing,MM.RankCheck} = nothing
    dependence::Type = MM.StaircaseDependence
end

function column_maxdegree(s::Solver, polys)
    if iszero(s.column_maxdegree)
        return sum(MP.maxdegree, polys) - length(polys) + 2
    else
        return s.column_maxdegree
    end
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
    # Type not concrete as different iteration might give different type if the user changes the option
    border::Union{Nothing,MM.BorderBasis}
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
    if !isnothing(s.border)
        println(io, s.border)
    end
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

_some_args(::Nothing) = tuple()
_some_args(arg) = (arg,)

function expand!(s::Iterator, sel::AbstractShiftsSelector)
    return expand!(s.matrix, sel; sparse_columns = s.solver.sparse_columns)
end

function expand!(s::Iterator, sel::FirstStandardNonSaturated)
    targets = MP.monomial_type(typeof(s.matrix))[]
    for mono in MM.standard_basis(s.border.dependence).monomials
        if length(targets) >= sel.max_num
            break
        end
        # If we don't check for `is_forever_trivial` then we might
        # think we are adding monomials but we add none and hence
        # we won't make any progress
        if !is_saturated(s.matrix, mono) &&
            !is_forever_trivial(s.matrix, mono)
            push!(targets, mono)
        end
    end
    if isempty(targets) && s.solver.print_level >= 1
        @info("No candidate to saturate")
        return 0
    end
    added = expand!(
        s.matrix,
        TargetColumns(targets);
        sparse_columns = s.solver.sparse_columns,
    )
    if added > 0 && s.solver.print_level >= 1
        @info("Added $added rows to saturate columns `$targets`")
    end
    return added
end

function expand!(s::Iterator, ::Type{ColumnDegrees})
    mindeg = maximum(MP.maxdegree, s.matrix.polynomials)
    maxdeg = column_maxdegree(s.solver, s.matrix.polynomials)
    deg = mindeg - 1
    added = 0
    while iszero(added) && deg < maxdeg
        deg += 1
        if deg == mindeg
            min = minimum(MP.maxdegree, s.matrix.polynomials)
            degs = min:deg
        else
            degs = deg:deg
        end
        added = expand!(s, ColumnDegrees(degs))
    end
    if added > 0 && s.solver.print_level >= 1
        @info("Added $added rows to complete columns up to degree $deg")
    end
    return added
end

function step!(s::Iterator, it)
    if s.solver.max_iter > 0 && size(s.stats, 1) >= s.solver.max_iter
        if s.solver.print_level >= 1
            @info("Reached iteration limit of $(s.solver.max_iter) iterations")
        end
        s.status = MOI.ITERATION_LIMIT
        return
    end
    added = expand!(s, it)
    if iszero(added)
        if s.solver.print_level >= 1
            @info("No row added")
        end
        s.status = MOI.OTHER_LIMIT
        return
    end
    if s.solver.print_level >= 3
        display(s)
    end
    null = LinearAlgebra.nullspace(s.matrix, _some_args(s.solver.rank_check)...)
    nullity = size(null.matrix, 2)
    try_solving = !s.solver.wait_for_gap || (!isempty(s.stats.nullity) && nullity == s.stats.nullity[end])
    if s.solver.trim_to_border || try_solving
        rank_check = something(s.solver.rank_check, MM.LeadingRelativeRankTol(1e-8))
        s.border = MM.BorderBasis{s.solver.dependence}(null, rank_check)
    end
    if try_solving
        sols = MM.solve(s.border)
        if !isnothing(sols)
            s.solutions = collect(sols)
            if isempty(s.solutions)
                if s.solver.print_level >= 1
                    @info("Infeasible system")
                end
                s.status = MOI.INFEASIBLE
            else
                if s.solver.print_level >= 1
                    @info("Found $(length(s.solutions)) real solution$(isone(length(s.solutions)) ? "" : "s")")
                end
                s.status = MOI.OPTIMAL
            end
        end
    end
    if s.solver.trim_to_border
        dep = s.border.dependence
        standard_and_border = MM.sub_basis(
            dep,
            findall(eachindex(dep.dependence)) do i
                if isnothing(MM._index(null.basis, dep.basis.monomials[i]))
                    return false
                end
                return MM.is_standard(dep.dependence[i]) || any(MP.variables(dep.basis)) do v
                    mono = dep.basis.monomials[i]
                    if MP.divides(v, mono)
                        q = MP.div_multiple(mono, v)
                        j = MM._index(dep.basis, q)
                        if !isnothing(j) && MM.is_standard(dep.dependence[j])
                            return true
                        end
                    end
                    return false
                end
            end
        )
        trimmed_null = null[standard_and_border]
        s.matrix = MacaulayMatrix(LinearAlgebra.nullspace(trimmed_null.matrix')' * standard_and_border.monomials)
        # To at least add the unshifted polynomials otherwise the next iteration
        # will do nothing
        expand!(s, it)
    end
    push!(s.stats, [nullity, size(s.matrix)...])
    return
end

include("saturation.jl")
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

export init, step!, solve!

include("cvp.jl")

include("H2/H2.jl")
export H2

end # module Macaulay
