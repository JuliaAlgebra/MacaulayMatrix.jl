module TestH2

using Test
using Macaulay

function test_Mauricio()
    system, vars = H2.composeSystemWalsh(H2.Mauricio(), 1)
    a1, b0, r0, r1 = vars
    @test system[1] == -b0 - r0 + 1
    @test system[2] == 2a1*r0 + a1 - 12b0 - r1 + 9
    @test system[3] == -a1^2 * r0 + 2a1 * r1 + 9a1 - 49b0 - 10
    @test system[4] == -a1^2 * r1 - 10a1 - 78b0
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

TestH2.runtests()
