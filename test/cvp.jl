module TestCVP

using LinearAlgebra
using Test
import Macaulay

function test_elim(T = Float64)
    M = T[
        1 0 1
        0 1 1
    ]
    x = Macaulay.eliminate_indices(M, [3])
    @test x / x[1] ≈ [1, -1, 0]
end

using TypedPolynomials
import MultivariatePolynomials as MP

function test_cvp(T = Float64)
    @polyvar x
    for α in [zero(T), one(T), T(2)]
        p = x^2 - 2x + α
        c, ε = Macaulay.cvp(p, 2)
        @test MP.maxdegree(c) == 1
        c /= MP.leading_coefficient(c)
        vars = MP.effective_variables(c)
        @test length(vars) == 1
        σ = first(vars)
        @test c ≈ one(T) * σ - α + one(T) atol=1e-3
        @test ε < 1e-3
    end
end

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$name", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
end

end

TestCVP.runtests()
