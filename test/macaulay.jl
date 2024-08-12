module TestMacaulayMatrix

using SparseArrays, Test
using TypedPolynomials
import MultivariateBases as MB
using MacaulayMatrix
using JuMP
import CSDP

function __test_monomial_ideal_generators(standard, generators)
    basis = MB.MonomialBasis(standard)
    output = @inferred monomial_ideal_generators(basis)
    @test output isa MB.MonomialBasis
    # FIXME equality between `MB.MonomialBasis` does not work
    # `MB.MonomialBasis(generators)` ensures the basis is sorted
    @test output.monomials == MB.MonomialBasis(generators).monomials
end

function _test_monomial_ideal_generators(generators)
    standard = filter(
        monomials(
            variables(generators),
            0:sum(maxdegree(generators, v) for v in variables(generators)),
        ),
    ) do mono
        return !any(generators) do g
            return divides(g, mono)
        end
    end
    __test_monomial_ideal_generators(standard, generators)
    return
end

function test_monomial_ideal_generators()
    @polyvar x y z
    _test_monomial_ideal_generators([x^0])
    _test_monomial_ideal_generators([x^3])
    _test_monomial_ideal_generators([x^3, x^2 * y^2, y^3])
    return _test_monomial_ideal_generators([
        x^7,
        x^3 * z,
        x^2 * z * y^2,
        y^8,
        y^3 * z^2,
        z^7,
    ])
end

# Taken from `macaulaylab.net/Tests/testMacaulayMatrix.m`
function test_MacaulayMatrix()
    @polyvar x y
    p = [1 + 2x + 3y + 4y^2, 5 + 6x]
    M = macaulay(p, 3)

    M_expected = [
        1 3 2 4 0 0 0 0 0 0
        5 0 6 0 0 0 0 0 0 0
        0 1 0 3 2 0 4 0 0 0
        0 5 0 0 6 0 0 0 0 0
        0 0 1 0 3 2 0 4 0 0
        0 0 5 0 0 6 0 0 0 0
        0 0 0 5 0 0 0 6 0 0
        0 0 0 0 5 0 0 0 6 0
        0 0 0 0 0 5 0 0 0 6
    ]

    @test sparse(M) == M_expected
end

# Taken from `macaulaylab.net/Database/Systems/dreesen1.m`
function dreesen1()
    @polyvar x y
    return [-x^2 + 2x * y + y^2 + 5x - 3y - 4.0, x^2 + 2x * y + y^2 - 1]
end

function _test_sols(sols, expected)
    @test length(sols) == length(expected)
    for exp in expected
        @test any(sol -> isapprox(sol, exp, rtol = 1e-3), sols)
    end
end

function test_dreesen1()
    ps = dreesen1()
    expected = [[4, -5], [3, -2], [0, -1], [1, 0]]
    @testset "d=$d" for d in 3:5
        @testset "solve_system $sparse_columns" for sparse_columns in
                                                    [false, true]
            @testset "trim_to_border $trim_to_border" for trim_to_border in
                                                          [false, true]
                solver = Solver(;
                    trim_to_border,
                    column_maxdegree = d,
                    sparse_columns,
                )
                s = MacaulayMatrix.SS.algebraic_set(ps, solver)
                sols = collect(s)
                _test_sols(sols, expected)
            end
        end
        @testset "psd_hankel" begin
            solver =
                optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)
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
    @test solve_system([q], column_maxdegree = 3, wait_for_gap = true) ===
          nothing
    _test_sols(solve_system([q], column_maxdegree = 3), [[exp]])
    for d in 4:4
        sols = solve_system([q], column_maxdegree = d)
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

TestMacaulayMatrix.runtests()
