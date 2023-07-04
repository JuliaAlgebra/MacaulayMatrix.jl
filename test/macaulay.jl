module TestMacaulay

using Test, TypedPolynomials
using Macaulay
using JuMP
import CSDP

# Taken from `macaulaylab.net/Tests/testmacaulay.m`
function test_macaulay()
    @polyvar x y
    p = [1 + 2x + 3y + 4y^2, 5 + 6x]
    M = macaulay(p, 3)

    M_expected = [
        5 0 6 0 0 0 0 0 0 0
        1 3 2 4 0 0 0 0 0 0
        0 5 0 0 6 0 0 0 0 0
        0 0 5 0 0 6 0 0 0 0
        0 1 0 3 2 0 4 0 0 0
        0 0 1 0 3 2 0 4 0 0
        0 0 0 5 0 0 0 6 0 0
        0 0 0 0 5 0 0 0 6 0
        0 0 0 0 0 5 0 0 0 6
    ]

    @test M == M_expected
end

# Taken from `macaulaylab.net/Database/Systems/dreesen1.m`
function dreesen1()
    @polyvar x y
    return [
        -x^2 + 2x * y + y^2 + 5x - 3y - 4,
         x^2 + 2x * y + y^2           - 1,
    ]
end

function _test_sols(sols, expected)
    @test length(sols) == length(expected)
    for exp in expected
        @test any(sol -> isapprox(sol, exp, rtol=1e-3), sols)
    end
end

function test_dreesen1() 
    ps = dreesen1()
    expected = [
        [4, -5],
        [3, -2],
        [0, -1],
        [1, 0],
    ]
    solver = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)
    @testset "d=$d" for d in 3:5
        @testset "solve_system" begin
            sols = solve_system(ps, d)
            _test_sols(sols, expected)
        end
        @testset "psd_hankel" begin
            sols = psd_hankel(ps, solver, d)
            if d == 3
                @test sols === nothing
            else
                _test_sols(sols, expected)
            end
        end
    end
end

function test_univariate()
    @polyvar x
    p = 3x^4 + 8x^3 - 6x^2 + 24x + 1
    q = differentiate(p, x)
    exp = -2.658967
    @test solve_system([q], 3) === nothing
    for d in 4:8
        sols = solve_system([q], d)
        _test_sols(sols, [[exp]])
    end
    solver = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)
    @test psd_hankel([q], solver, 3) === nothing
    @testset "d=$d" for d in 4:8
        sols = psd_hankel([q], solver, d)
        if isodd(d)
            @test sols === nothing # FIXME
        else
            _test_sols(sols, [[exp]])
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

TestMacaulay.runtests()
