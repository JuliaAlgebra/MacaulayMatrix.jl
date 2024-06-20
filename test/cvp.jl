module TestCVP

using LinearAlgebra
using Test
import MacaulayMatrix
using Random

function test_elim(T = Float64)
    Random.seed!(0)
    M = T[
        1 0 1
        0 1 1
    ]
    x = MacaulayMatrix.eliminate_indices(M, [3])
    @test x / x[1] ≈ [1, -1, 0] rtol = 1e-6
    x = MacaulayMatrix.eliminate_indices(M, [3], 2)
    @test x / x[1] ≈ [1, -1, 0] rtol = 1e-6
    x = MacaulayMatrix.eliminate_indices(M, [3], 2, solver = MacaulayMatrix.Manopt.LevenbergMarquardt)
    @test x / x[1] ≈ [1, -1, 0] rtol = 1e-6
    M = T[
        1 0 1 0
        0 1 0 1
        0 0 1 -1
    ]
    x = MacaulayMatrix.eliminate_indices(M, [3, 4], 2; solver = MacaulayMatrix.Manopt.LevenbergMarquardt)
    @test x / x[1] ≈ [1, -1, 0, 0] rtol = 1e-6
end

using TypedPolynomials
import MultivariatePolynomials as MP

function test_cvp(T = Float64)
    Random.seed!(0)
    @polyvar x
    @testset "$ellipsoid" for ellipsoid in [true, false]
        solvers = Any[]
        if ellipsoid
            # Failing on ci for 32-bits with non-ellipsoid
            push!(solvers, MacaulayMatrix.Manopt.GradientDescentState)
        else
            # Missing https://github.com/JuliaManifolds/Manifolds.jl/blob/106a8f35f25c5e0501b7b9ec7cf6456e4effe947/src/manifolds/Sphere.jl#L254-L265
            # for ellispoid to work with LevenbergMarquardt
            push!(solvers, MacaulayMatrix.Manopt.LevenbergMarquardt)
        end
        rank_checks = Any[nothing]
        if !ellipsoid
            push!(rank_checks, MacaulayMatrix.MM.LeadingRelativeRankTol(Base.rtoldefault(Float64)))
        end
        @testset "$(string(solver)[1:8])" for solver in solvers
            @testset "$(string(rank_check)[1:4])" for rank_check in rank_checks
                @testset "$α" for α in [zero(T), one(T), T(2)]
                    p = x^2 - 2x + α
                    c, ε = MacaulayMatrix.cvp(p, 2; ellipsoid, solver, rank_check)
                    @test MP.maxdegree(c) == 1
                    c /= MP.leading_coefficient(c)
                    vars = MP.effective_variables(c)
                    @test length(vars) == 1
                    σ = first(vars)
                    @test c ≈ one(T) * σ - α + one(T) atol=1e-3
                    @test ε < 1e-3
                end
            end
        end
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
