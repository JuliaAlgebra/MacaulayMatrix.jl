module Macaulay

import LinearAlgebra, SparseArrays
import MultivariatePolynomials as MP
import JuMP
import MultivariateMoments

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
    nM, cM, U = MultivariateMoments.lowrankchol(
        H,
        MultivariateMoments.SVDChol(),
        MultivariateMoments.LeadingRelativeRankTol(1e-6),
    )
    return realization_observability(U')
end

function realization_hankel(H::LinearAlgebra.Symmetric, monos)
    M = MultivariateMoments.MomentMatrix(H, monos)
    η = MultivariateMoments.extractatoms(M, 1e-4)
    if isnothing(η)
        return
    else
        return [η.atoms[i].center for i in eachindex(η.atoms)]
    end
end

function moment_matrix(Z::AbstractMatrix, solver, d, monos)
    model = JuMP.Model(solver)
    JuMP.@variable(model, b[1:size(Z, 2)])
    JuMP.@constraint(model, sum(b) == 1)
    Zb = Z * b
    gram_monos = MP.monomials(MP.variables(monos), 0:d)
    H = LinearAlgebra.Symmetric([
        Zb[findfirst(isequal(gram_monos[i] * gram_monos[j]), monos)]
        for i in eachindex(gram_monos), j in eachindex(gram_monos)
    ])
    JuMP.@constraint(model, H in JuMP.PSDCone())
    JuMP.optimize!(model)
    if JuMP.termination_status(model) == JuMP.MOI.INFEASIBLE
        return
    elseif JuMP.termination_status(model) != JuMP.MOI.OPTIMAL
        error(string(JuMP.solution_summary(model)))
    end
    return LinearAlgebra.Symmetric(JuMP.value.(H)), gram_monos
end

function psd_hankel(Z::AbstractMatrix, solver, d, monos)
    H = moment_matrix(Z, solver, d, monos)
    if H === nothing
        return nothing
    end
    return realization_hankel(H...)
end

function psd_hankel(polynomials::AbstractVector{<:MP.AbstractPolynomialLike{T}}, solver, maxdegree) where {T}
    M, monos = macaulay_monomials(polynomials, maxdegree)
    Z = LinearAlgebra.nullspace(Matrix(M))
    return psd_hankel(Z, solver, div(maxdegree, 2), monos)
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
            lm = MP.leadingmonomial(p)
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

# TODO implement sieve
function standard_monomials(Z, tol = 1e-10)
    list = Int[]
    old_rank = 0
    for k in axes(Z, 1)
        new_rank = LinearAlgebra.rank(Z[1:k, :], tol)
        if new_rank > old_rank
            push!(list, k)
        end
        old_rank = new_rank
        if new_rank == size(Z, 2)
            break
        end
    end
    return list
end

num_monomials(d, n) = binomial(n + d, d)

function gap_zone_standard_monomials(c, maxdegree, n)
    ma = 0
    gapsize = 0
    dgap = nothing
    for d in 0:maxdegree
        num = num_monomials(d, n)
        r = count(c .<= num) # FIXME findsortedfirst
        if ma == r
            gapsize += 1;
            dgap = d - gapsize + 1;
        else
            if isnothing(dgap)
                ma = r
            else
                break
            end
        end
    end     
    return dgap, ma, gapsize
end

function shiftnullspace(Z, dgap, srows, monos)
    S = Z[srows, :]
    Sx = [Z[[findfirst(isequal(monos[row]* x), monos) for row in srows], :] for x in MP.variables(monos)]
    pS = LinearAlgebra.pinv(S)
    Sx = [pS * S for S in Sx]
    return MultivariateMoments.SemialgebraicSets.solvemultiplicationmatrices(
        Sx,
        MultivariateMoments.SemialgebraicSets.ReorderedSchurMultiplicationMatricesSolver{Float64}(),
    )

    return MultivariateMoments.solve
    eig = [LinearAlgebra.eigen(S).values for S in Sx]
    return [[eig[i][j] for i in eachindex(eig)] for j in eachindex(eig[1])]
end

# Inspired from `macaulaylab.net/Code/solvesystemnullspace.m`
function solve_system(polynomials::AbstractVector{<:MP.AbstractPolynomialLike{T}}, maxdegree) where {T}
    mindegree = maximum(MP.maxdegree, polynomials)
    n = MP.nvariables(polynomials)
    nullities = zeros(Int, maxdegree)
    ma = mb = nothing
    Z = nothing
    dgap = nothing
    srows = nothing
    monos = nothing
    for d in mindegree:maxdegree
        M, monos = macaulay_monomials(polynomials, d)
        Z = LinearAlgebra.nullspace(Matrix(M))
        nullities[d] = size(Z, 2)
        if d > mindegree && nullities[d] == nullities[d - 1]
            sols = realization_observability(Z, monos)
            if !isnothing(sols)
                return sols
            end
        end
    end
    return
end

function realization_observability(Z, monos)
    n = MP.nvariables(monos)
    d = MP.maxdegree(monos)
    srows = standard_monomials(Z)
    dgap, ma, gapsize = gap_zone_standard_monomials(srows, d, n)
    srows = srows[1:ma]
    if gapsize < 1
        return
    end
    mb = size(Z, 2)
    if mb == ma 
        W = Z
    else
        error("Column compression not supported yet")
    end

    # Solve the system:
    return shiftnullspace(W, dgap, srows, monos)

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

end # module Macaulay
