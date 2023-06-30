using Test #src
using DynamicPolynomials
@polyvar x y
p = x^3 - x^2 + 2x*y -y^2 + y^3
using SumOfSquares
S = @set x >= 0 && y >= 0 && x + y >= 1
p(x=>1, y=>0), p(x=>1//2, y=>1//2), p(x=>0, y=>1)


import Macaulay
model = Model(PolyJuMP.KKT.Optimizer)
set_attribute(model, "algebraic_solver", Macaulay.Solver(3))
@variable(model, a >= 0)
@variable(model, b >= 0)
@constraint(model, a + b >= 1)
@NLobjective(model, Min, a^3 - a^2 + 2a*b - b^2 + b^3)
optimize!(model)

# As we can see below, the termination status is `LOCALLY_SOLVED` and not of `OPTIMAL`
# because Ipopt only guarantees **local** optimality.

solution_summary(model)

for i in 1:result_count(model)
    println("$i: $(value(a, result = i)), $(value(b, result = i))")
end

# Indeed, the solution found is not globally optimal:

@test value(a) ≈ 0.5 rtol=1e-5 #src
@test value(b) ≈ 0.5 rtol=1e-5 #src
value(a), value(b)

# Note that the problem can be written equivalently as follows using [registered functions](https://jump.dev/JuMP.jl/stable/manual/nlp/#Register-a-function).
# The difference is that the gradient and hessian will be computed via the *Symbolic Differentiation* provided
# by MultivariatePolynomials instead of JuMP's *Automatic Differentiation*:

f(a, b) = p(x => a, y => b)
∇p = differentiate(p, [x, y])
function ∇f(g, a, b)
    for i in eachindex(g)
        g[i] = ∇p[i](x => a, y => b)
    end
end
∇²p = differentiate(∇p, [x, y])
function ∇²f(H, a, b)
    for j in axes(∇²p, 2)
        for i in j:size(∇²p, 1)
            H[i, j] = ∇²p[i, j](x => a, y => b)
        end
    end
end
using Ipopt
gmodel = Model(Ipopt.Optimizer)
@variable(gmodel, a >= 0)
@variable(gmodel, b >= 0)
@constraint(gmodel, a + b >= 1)
register(gmodel, :f, 2, f, ∇f, ∇²f)
@NLobjective(gmodel, Min, f(a, b))
optimize!(gmodel)

# Even if we have the algebraic expressions of gradient and hessian,
# Ipopt is not using these symbolic expressions but only local information
# hence it can still only provide local guarantees:

@test termination_status(gmodel) == MOI.LOCALLY_SOLVED #src
@test objective_value(gmodel) ≈ 0.25 rtol=1e-5 #src
solution_summary(gmodel)

#
