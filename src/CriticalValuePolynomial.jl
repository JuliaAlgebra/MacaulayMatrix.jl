module CriticalValuePolynomial

import MultivariatePolynomials as MP
import SparseArrays

export macaulay

function macaulay(polynomials::AbstractVector{<:MP.AbstractPolynomialLike{T}}, maxdegree) where {T}
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
    return SparseArrays.sparse(I, J, K, row, length(monos))
end

end # module CriticalValuePolynomial
