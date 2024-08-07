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

function MM.moment_matrix(null::MM.MacaulayNullspace, solver, d; print_level = 1, T = Float64)
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
            return convert(JuMP.GenericAffExpr{T,JuMP.GenericVariableRef{T}}, b[r+inf_idx])
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
    return MM.moment_matrix(Z, solver, div(maxdegree, 2); kws...)
end

function psd_hankel(args...)
    H = MM.moment_matrix(args...)
    if H === nothing
        return nothing
    end
    return realization_hankel(H)
end
