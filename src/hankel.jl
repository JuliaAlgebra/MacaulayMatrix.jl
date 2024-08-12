export solutions,
    errors, support_error, cheat_rank, cheat_nullspace, cheat_system

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

function MM.moment_matrix(
    null::MM.MacaulayNullspace,
    solver,
    d = div(MP.maxdegree(null.basis), 2);
    print_level = 1,
    T = Float64,
)
    # TODO Newton polytope
    vars = MP.variables(null.basis.monomials)
    monos = MP.monomials(vars, 0:2d)
    model = JuMP.GenericModel{T}(solver)
    r = size(null.matrix, 2)
    # Number of roots at infinity that were left out by not adding these as columns
    null = null[[mono for mono in monos if mono in null.basis.monomials]]
    num_inf = length(monos) - size(null.matrix, 1)
    JuMP.@variable(model, b[1:(r+num_inf)])
    JuMP.@constraint(model, sum(b) == 1)
    inf_idx = 0
    Zb = map(monos) do mono
        idx = MM._monomial_index(null.basis.monomials, mono)
        if isnothing(idx)
            inf_idx += 1
            return convert(
                JuMP.GenericAffExpr{T,JuMP.GenericVariableRef{T}},
                b[r+inf_idx],
            )
        else
            return null.matrix[idx, :]' * b[1:r]
        end
    end
    @assert inf_idx == num_inf
    gram_monos = MP.monomials(vars, 0:d)
    H = LinearAlgebra.Symmetric([
        Zb[findfirst(isequal(gram_monos[i] * gram_monos[j]), monos)] for
        i in eachindex(gram_monos), j in eachindex(gram_monos)
    ])
    JuMP.@constraint(model, H in JuMP.PSDCone())
    JuMP.optimize!(model)
    if print_level >= 1
        @info(
            "Terminated with $(JuMP.termination_status(model)) ($(JuMP.raw_status(model))) in $(JuMP.solve_time(model)) seconds."
        )
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

function MM.moment_matrix(
    polynomials::AbstractVector{<:MP.AbstractPolynomialLike{T}},
    solver,
    maxdegree;
    kws...,
) where {T}
    Z = LinearAlgebra.nullspace(macaulay(polynomials, maxdegree))
    return MM.moment_matrix(Z, solver; kws...)
end

function psd_hankel(args...)
    H = MM.moment_matrix(args...)
    if H === nothing
        return nothing
    end
    return realization_hankel(H)
end

function solutions(M::MM.MomentMatrix, s = MM.ShiftNullspace())
    sols = []
    for r in 1:length(M.basis)
        η = MM.atomic_measure(M, MM.FixedRank(r), s)
        if !isnothing(η)
            for atom in η.atoms
                push!(sols, r => atom)
            end
        end
    end
    return sols
end

function sol_diracs(ν, vars, sols)
    return [MM.dirac(ν.basis.monomials, vars => sol) for sol in sols]
end

function errors(ν::MM.MomentMatrix, vars, sols)
    diracs = sol_diracs(ν, vars, sols)
    S = LinearAlgebra.svd(MM.value_matrix(ν))
    errors = fill(NaN, size(S.U, 2))
    for i in size(S.U, 2):-1:1
        u = S.U[:, i]
        errors[i] =
            maximum(diracs, init = i == size(S.U, 2) ? 0 : errors[i+1]) do dirac
                return abs(LinearAlgebra.dot(dirac.a, u))
            end
    end
    return errors
end

function support_error(ν::MM.MomentMatrix, vars, sols)
    return maximum(SS.equalities(ν.support), init = 0) do eq
        return maximum(sols) do sol
            return abs(eq(vars => sol))
        end
    end
end

function cheat_rank(ν::MM.MomentMatrix, vars, sols, rank_check)
    return MM.rank_from_singular_values(errors(ν, vars, sols), rank_check)
end

"""
    cheat_nullspace(ν::MM.MomentMatrix, args...)

Return the nullspace of `ν` as matrix for each each row
is an element of the nullspace.
"""
function cheat_nullspace(ν::MM.MomentMatrix, args...)
    r = cheat_rank(ν, args...)
    return LinearAlgebra.nullspace(ν, MM.FixedRank(r))
end

"""
    cheat_system(ν::MM.MomentMatrix, args...)

Return the nullspace of `ν` as an vector of polynomials.
"""
function cheat_system(ν::MM.MomentMatrix, args...)
    return cheat_nullspace(ν, args...).polynomials
end

function LinearAlgebra.nullspace(ν::MM.MomentMatrix, rank_check::MM.RankCheck)
    S = LinearAlgebra.svd(MM.value_matrix(ν))
    r = MM.rank_from_singular_values(S.S, rank_check)
    return LazyMatrix(S.U[:, (r+1):end]' * ν.basis.monomials)
end

function LinearAlgebra.nullspace(ν::MM.MomentMatrix, rank_check::MM.RankCheck, solver::MM.ShiftNullspace; kws...)
    null = MM.MacaulayNullspace(ν, rank_check)
    border = MM.BorderBasis{MM.StaircaseDependence}(null, solver.check)
    std = MM.standard_basis(border.dependence; trivial = false)
    dep = MM.dependent_basis(border.dependence)
    ν.support = SS.algebraic_set(dep.monomials - border.matrix' * std.monomials)
    return LinearAlgebra.nullspace(ν; kws...)
end

function LinearAlgebra.nullspace(ν::MM.MomentMatrix, rank_check::MM.RankCheck, solver::MM.Echelon; kws...)
    MM.compute_support!(ν, rank_check, solver)
    return LinearAlgebra.nullspace(ν; kws...)
end

function clean(α::Real; tol = 1e-8)
    if abs(α) < tol
        return zero(α)
    else
        r = round(α)
        if abs(α - r) < tol
            return r
        end
        return α
    end
end

function clean(p::MP.AbstractPolynomialLike; kws...)
    return MP.map_coefficients(α -> clean(α; kws...), p)
end

function LinearAlgebra.nullspace(ν::MM.MomentMatrix; kws...)
    return LazyMatrix(clean.(SS.equalities(ν.support); kws...))
end
