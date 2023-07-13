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
