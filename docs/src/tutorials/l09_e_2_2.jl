using Test #src
using DynamicPolynomials
@polyvar x y
p = x^3 - x^2 + 2x*y -y^2 + y^3
using SumOfSquares
S = @set x >= 0 && y >= 0 && x + y >= 1
p(x=>1, y=>0), p(x=>1//2, y=>1//2), p(x=>0, y=>1)

@polyvar x[1:2]
@polyvar σ[1:3]
S = @set x[1] * σ[1] == 0 &&
         x[2] * σ[2] == 0 &&
         x[2] * σ[3] + x[1] * σ[3] == σ[3] &&
         2.0x[2] - 2x[1] - σ[3]^2 - σ[1]^2 + 3x[1]^2 == 0 &&
        -2x[2] + 2x[1] - σ[3]^2 - σ[2]^2 + 3x[2]^2 == 0 solver = Macaulay.Solver(6, 3)
display(S)
S = @set x[1] * σ[1] == 0 &&
         x[2] * σ[2] == 0 &&
         x[1] + x[2] == 1 &&
         x[2] * σ[3] + x[1] * σ[3] == σ[3] &&
         2.0x[2] - 2x[1] - σ[3]^2 - σ[1]^2 + 3x[1]^2 == 0 &&
        -2x[2] + 2x[1] - σ[3]^2 - σ[2]^2 + 3x[2]^2 == 0 solver = Macaulay.Solver(6, 3)
collect(S)

S = @set x[1] * σ[1] == 0 &&
         x[2] * σ[2] == 0 &&
         x[1] + x[2] == 1 &&
         2.0x[2] - 2x[1] - σ[1]^2 + 3x[1]^2 == 1 &&
        -2x[2] + 2x[1] - σ[2]^2 + 3x[2]^2 == 1 solver = Macaulay.Solver(6, 3)
collect(S)


import Macaulay
model = Model(PolyJuMP.KKT.Optimizer)
@variable(model, a >= 0)
@variable(model, b >= 0)
@constraint(model, a + b >= 1)
@NLobjective(model, Min, a^3 - a^2 + 2a*b - b^2 + b^3)
set_attribute(model, "algebraic_solver", Macaulay.Solver(6, 3))
optimize!(model)

# As we can see below, the termination status is `LOCALLY_SOLVED` and not of `OPTIMAL`
# because Ipopt only guarantees **local** optimality.

solution_summary(model)

for i in 1:result_count(model)
    println("$i: $(value(a, result = i)), $(value(b, result = i))")
end
