export solutions, errors, support_error
export cheat_rank, cheat_nullspace, cheat_system
export Hankel, Explicit

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

function _check_status(model, status)
    if JuMP.MOI.get(model, status) != JuMP.MOI.FEASIBLE_POINT ||
       JuMP.termination_status(model) != JuMP.MOI.OPTIMAL
        message = string(JuMP.solution_summary(model))
        if JuMP.MOI.get(model, status) == JuMP.MOI.NO_SOLUTION
            error(message)
        else
            @warn(message)
        end
    end
end

function MM.moment_matrix(
    null::MM.MacaulayNullspace,
    solver,
    d = div(MP.maxdegree(null.basis), 2);
    print_level = 1,
    T = Float64,
)
    if !MP.isconstant(null.basis.monomials[1])
        error("The constant monomial is not part of the basis.")
    end
    # TODO Newton polytope
    vars = MP.variables(null.basis.monomials)
    monos = MP.monomials(vars, 0:2d)
    model = JuMP.GenericModel{T}(solver)
    r = size(null.matrix, 2)
    # Number of roots at infinity that were left out by not adding these as columns
    null = null[[mono for mono in monos if mono in null.basis.monomials]]
    num_inf = length(monos) - size(null.matrix, 1)
    JuMP.@variable(model, b[1:(r+num_inf)])
    inf_idx = 0
    missing_monos = eltype(null.basis.monomials)[]
    Zb = map(monos) do mono
        idx = MM._monomial_index(null.basis.monomials, mono)
        if isnothing(idx)
            inf_idx += 1
            push!(missing_monos, mono)
            return convert(
                JuMP.GenericAffExpr{T,JuMP.GenericVariableRef{T}},
                b[r+inf_idx],
            )
        else
            return null.matrix[idx, :]' * b[1:r]
        end
    end
    if !isempty(missing_monos)
        @warn("Missing monomials $(missing_monos) in Macaulay nullspace, leaving them free to take any value in the moment matrix then.")
    end
    JuMP.@constraint(model, Zb[1] == 1)
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
    else
        _check_status(model, JuMP.MOI.PrimalStatus())
    end
    H = LinearAlgebra.Symmetric(JuMP.value.(H))
    return MM.MomentMatrix(H, gram_monos)
end

struct Hankel end
struct Explicit end

# Moment problem should be a pure feasibility problem
# to help maximize the rank.
# min γ
# γ + ∑ λ_i(x) p_i(x) is SOS
#
# ⟨μ, p_i(x) λ(x)⟩ = 0 ∀i, λ (Localization matrix)
# ⟨μ, 1⟩ = 1
function MM.moment_matrix(
    polynomials::AbstractVector{<:MP.AbstractPolynomialLike},
    solver,
    maxdegree,
    ::Hankel;
    kws...,
)
    Z = LinearAlgebra.nullspace(macaulay(polynomials, maxdegree))
    return MM.moment_matrix(Z, solver; kws...)
end

function MM.moment_matrix(
    polynomials::AbstractVector,
    solver,
    maxdegree::Integer;
    kws...,
)
    return MM.moment_matrix(polynomials, solver, maxdegree, Hankel(); kws...)
end

function MM.moment_matrix(
    polynomials::AbstractVector,
    solver,
    maxdegree::Integer,
    ::Explicit;
    kws...,
)
    model = JuMP.Model(solver)
    JuMP.@variable(model, γ)
    JuMP.@objective(model, Max, γ)
    vars = MP.variables(polynomials)
    JuMP.@variable(
        model,
        λ[i in eachindex(polynomials)],
        SumOfSquares.PolyJuMP.Poly(
            MP.monomials(vars, 0:(maxdegree-MP.maxdegree(polynomials[i]))),
        ),
    )
    con_ref = JuMP.@constraint(
        model,
        LinearAlgebra.dot(polynomials, λ) - γ in SumOfSquares.SOSCone()
    )
    JuMP.optimize!(model)
    if JuMP.termination_status(model) == JuMP.MOI.DUAL_INFEASIBLE
        return
    else
        _check_status(model, JuMP.MOI.DualStatus())
    end
    return MM.moment_matrix(con_ref)
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

function LinearAlgebra.nullspace(ν::MM.MomentMatrix, rank_check)
    S = LinearAlgebra.svd(MM.value_matrix(ν))
    r = MM.rank_from_singular_values(S.S, default_rank_check(rank_check))
    return LazyMatrix(S.U[:, (r+1):end]' * ν.basis.monomials)
end

function LinearAlgebra.nullspace(
    ν::MM.MomentMatrix,
    solver::MM.ShiftNullspace,
    rank_check = nothing;
    kws...,
)
    null = MM.image_space(ν, default_rank_check(rank_check))
    border = MM.BorderBasis{MM.StaircaseDependence}(null, solver.check)
    std = MM.standard_basis(border.dependence; trivial = false)
    dep = MM.dependent_basis(border.dependence)
    ν.support = SS.algebraic_set(dep.monomials - border.matrix' * std.monomials)
    return LinearAlgebra.nullspace(ν; kws...)
end

function LinearAlgebra.nullspace(
    ν::MM.MomentMatrix,
    solver::MM.Echelon,
    rank_check = nothing;
    kws...,
)
    MM.compute_support!(ν, default_rank_check(rank_check), solver)
    return LinearAlgebra.nullspace(ν; kws...)
end

function LinearAlgebra.nullspace(ν::MM.MomentMatrix)
    return LazyMatrix(SS.equalities(ν.support))
end
