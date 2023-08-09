"""
    @enum(ShiftStatus, NOT_INCLUDED, INCLUDED, NOT_REDUNDANT, REDUNDANT)

* `NOT_INCLUDED` : Not included yet in the rows of the Macaulay matrix
* `INCLUDED` : Already included in the rows of the Macaulay matrix but redundancy unknown
* `NOT_REDUNDANT` : Already included in the rows of the Macaulay matrix and no redundancy with the other rows included
* `REDUNDANT` : Redundant with another row already included in the Macaulay matrix
"""
@enum(ShiftStatus, NOT_INCLUDED, INCLUDED, NOT_REDUNDANT, REDUNDANT)

mutable struct MacaulayMatrix{
    T,
    P<:MP.AbstractPolynomialLike,
    V<:AbstractVector{P},
    B,
}
    polynomials::V
    row_shifts::B
    shift_statuses::Vector{Vector{ShiftStatus}}
    column_basis::B
end

function MacaulayMatrix(
    polynomials::V,
) where {
    T,
    P<:MP.AbstractPolynomialLike{T},
    V<:AbstractVector{P},
}
    M = MP.monomial_type(polynomials)
    B = MB.MonomialBasis{M,MP.monomial_vector_type(M)}
    return MacaulayMatrix{T,P,V,B}(
        polynomials,
        MB.empty_basis(B),
        Vector{ShiftStatus}[],
        MB.empty_basis(B),
    )
end

function Base.show(io::IO, M::MacaulayMatrix)
    num_rows, num_cols = size(M)
    println(io, "$(num_rows)×$(num_cols) Macaulay matrix for polynomials:")
    for poly in M.polynomials
        println(io, "  ", poly)
    end
    println(io, "The row shifts are:")
    println(io, M.row_shifts)
    println(io, "The column basis is:")
    println(io, M.column_basis)
end

function Base.size(M::MacaulayMatrix, i::Int)
    if i == 1
        return sum(M.shift_statuses, init = 0) do statuses
            return count(statuses) do s
                return s == INCLUDED || s == NOT_REDUNDANT
            end
        end
    elseif i == 2
        return length(M.column_basis)
    else
        return 1
    end
end

function Base.size(M::MacaulayMatrix)
    return (size(M, 1), size(M, 2))
end

function _merge_bases(a::MB.MonomialBasis, b::MB.MonomialBasis)
    return MB.MonomialBasis(MP.merge_monomial_vectors([
        a.monomials,
        b.monomials,
    ]))
end

abstract type AbstractShiftsSelector end

struct ColumnDegrees{V} <: AbstractShiftsSelector
    degrees::V
end

function select_shifts(vars, poly, d::ColumnDegrees)
    deg = MP.maxdegree(poly)
    if deg > maximum(d.degrees)
        return nothing
    else
        return MP.monomials(vars, max.(0, d.degrees .- deg))
    end

end

struct TargetColumns{M} <: AbstractShiftsSelector
    targets::Vector{M}
end

function select_shifts(vars, poly, d::TargetColumns{M}) where {M}
    mono = MP.leading_monomial(poly)
    return M[
        MP.div_multiple(target, mono)
        for target in d.targets
        if MP.divides(mono, target)
    ]
end

struct FixedShifts{M} <: AbstractShiftsSelector
    shifts::Vector{M}
end

function select_shifts(vars, poly, d::FixedShifts)
    return d.shifts
end

"""
    expand!(M::MacaulayMatrix, degs; sparse_columns::Bool = true)

Add the row shifts for which the maxdegree of the shifted polynomial would
be in `degs`. If `sparse_columns` is `false` then all monomials of degree in
`degs` will be added to the columns. Otherwise, only the columns that
correspond to the monomial of one of the shifted polynomials will be added.
"""
function expand!(M::MacaulayMatrix, shifts_selector; sparse_columns::Bool = true)
    vars = MP.variables(M.polynomials)
    MT = MP.monomial_type(eltype(M.polynomials))
    row_monos_to_add = Dict{MT,Vector{ShiftStatus}}()
    if sparse_columns
        col_monos_to_add = Set{MT}()
    else
        col_maxdeg = MP.maxdegree(M.column_basis.monomials)
    end
    added = 0
    for (j, poly) in enumerate(M.polynomials)
        selected_shifts = select_shifts(vars, poly, shifts_selector)
        if isnothing(selected_shifts)
            continue
        end
        for shift in selected_shifts
            i = MM._monomial_index(M.row_shifts.monomials, shift)
            if isnothing(i)
                if !haskey(row_monos_to_add, shift)
                    row_monos_to_add[shift] = fill(NOT_INCLUDED, length(M.polynomials))
                end
                statuses = row_monos_to_add[shift]
            else
                statuses = M.shift_statuses[i]
            end
            if statuses[j] == NOT_INCLUDED
                added += 1
                statuses[j] = INCLUDED
                for mono in MP.monomials(poly)
                    col = shift * mono
                    if isnothing(MM._monomial_index(M.column_basis.monomials, col))
                        if sparse_columns
                            if !(col in col_monos_to_add)
                                push!(col_monos_to_add, col)
                            end
                        else
                            col_maxdeg = max(col_maxdeg, MP.degree(col))
                        end
                    end
                end
            end
        end
    end
    if sparse_columns
        if !isempty(col_monos_to_add)
            M.column_basis = _merge_bases(
                M.column_basis,
                MB.MonomialBasis(collect(col_monos_to_add)),
            )
        end
    else
        M.column_basis = MB.MonomialBasis(MP.monomials(vars, 0:col_maxdeg))
    end
    if !isempty(row_monos_to_add)
        old_shifts = M.row_shifts
        M.row_shifts = _merge_bases(
            M.row_shifts,
            MB.MonomialBasis(collect(keys(row_monos_to_add))),
        )
        old_statuses = M.shift_statuses
        M.shift_statuses = Vector{Vector{ShiftStatus}}(undef, length(M.row_shifts))
        for (shift, statuses) in zip(old_shifts.monomials, old_statuses)
            i = MM._monomial_index(M.row_shifts.monomials, shift)
            M.shift_statuses[i] = statuses
        end
        for (shift, statuses) in row_monos_to_add
            i = MM._monomial_index(M.row_shifts.monomials, shift)
            M.shift_statuses[i] = statuses
        end
    end
    return added
end

function macaulay(polynomials, maxdegree; kws...)
    M = MacaulayMatrix(polynomials)
    expand!(M, ColumnDegrees(0:maxdegree); kws...)
    return M
end

function SparseArrays.sparse(M::MacaulayMatrix{T}) where {T}
    column = Dict(M.column_basis.monomials[i] => i for i in eachindex(M.column_basis.monomials))
    row = 0
    I = Int[]
    J = Int[]
    K = T[]
    for (shift, statuses) in zip(M.row_shifts.monomials, M.shift_statuses)
        for (j, status) in enumerate(statuses)
            if status == INCLUDED || status == NOT_REDUNDANT
                row += 1
                for t in MP.terms(M.polynomials[j])
                    push!(I, row)
                    push!(J, column[shift * MP.monomial(t)])
                    push!(K, MP.coefficient(t))
                end
            end
        end
    end
    @assert row == size(M, 1)
    return SparseArrays.sparse(I, J, K, row, size(M, 2))
end

function _nullspace(
    M::Matrix,
    # This corresponds to the default of `LinearAlgebra.nullspace`
    rank_check=MM.LeadingRelativeRankTol(min(size(M)...) * eps(real(float(oneunit(eltype(M)))))),
)
    m, n = size(M)
    if iszero(m) || iszero(n)
        T = LinearAlgebra.eigtype(eltype(M))
        Z = Matrix{T}(LinearAlgebra.I, n, n)
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

function LinearAlgebra.nullspace(M::MacaulayMatrix, args...)
    Δt = @elapsed begin
        S = SparseArrays.sparse(M)
        Z, accuracy = _nullspace(S, args...)
    end
    @info("Nullspace of dimensions $(size(Z)) computed from Macaulay matrix of dimension $(size(S)) in $Δt seconds.")
    return MM.MacaulayNullspace(Z, M.column_basis, accuracy)
end
