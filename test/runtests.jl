module TestMacaulay

using Test, TypedPolynomials
using CriticalValuePolynomial

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
