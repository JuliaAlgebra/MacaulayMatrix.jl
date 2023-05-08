module H2

import MutableArithmetics as MA
using DynamicPolynomials
using LinearAlgebra
import QuadGK
import PolynomialRoots
import TypedTables
import ControlSystemsBase
import Polynomials

abstract type SISOTransferFunction end

struct ContinuousSISOTransferFunction{T} <: SISOTransferFunction
    num::Vector{T}
    den::Vector{T}
end

struct DiscreteSISOTransferFunction{T} <: SISOTransferFunction
    num::Vector{T}
    den::Vector{T}
end

"""
Second order H2 optimal approximation of linear dynamical systems,
M.I. Ahmad, M. Frangos, I.M. Jaimoukha
(2011) - Proc. of the 18th IFAC World Congress, Milano, Italy.
"""
F1() = ContinuousSISOTransferFunction(
    [-2.9239, -39.5525, -97.5270, -147.1508],
    [1, 11.9584, 43.9119, 73.6759, 44.3821],
)

"""
Second order H2 optimal approximation of linear dynamical systems,
M.I. Ahmad, M. Frangos, I.M. Jaimoukha
(2011) - Proc. of the 18th IFAC World Congress, Milano, Italy.
"""
F2() = ContinuousSISOTransferFunction(
    [-1.2805, -6.2266, -12.8095, -9.3373],
    [1, 3.1855, 8.9263, 12.2936, 3.1987],
)

"""
Second order H2 optimal approximation of linear dynamical systems,
M.I. Ahmad, M. Frangos, I.M. Jaimoukha
(2011) - Proc. of the 18th IFAC World Congress, Milano, Italy.
"""
F3() = ContinuousSISOTransferFunction(
    [-1.3369, -4.8341, -47.5819, -42.7285],
    [1, 17.0728, 84.9908, 122.4400, 59.9309],
)

"""
Second order H2 optimal approximation of linear dynamical systems,
M.I. Ahmad, M. Frangos, I.M. Jaimoukha
(2011) - Proc. of the 18th IFAC World Congress, Milano, Italy.
"""
F6() = ContinuousSISOTransferFunction(
    [2, 11.5, 57.75, 178.625, 345.5, 323.625, 94.5],
    [1, 10, 46, 130, 239, 280, 194, 60.0],
)

"""
A New Algorithm for L_2 Optimal Model Reduction,
J.T. Spanos, M.H. Milman, D.L. Mingori
(1992) - Automatica.
"""
Spanos_1992_1() = ContinuousSISOTransferFunction(
    [1, 15, 50],
    [1, 5, 33, 79, 50],
)

"""
A New Algorithm for L_2 Optimal Model Reduction,
J.T. Spanos, M.H. Milman, D.L. Mingori
(1992) - Automatica.
"""
Spanos_1992_3() = ContinuousSISOTransferFunction(
    [0, 0, 0, 0.00001, 0.0110, 1],
    [1, 0.2220, 22.1242, 3.5445, 122.4433, 11.3231, 11.1100],
)


"""
H2 Model Reduction for Large-Scale Linear Dynamical Systems,
S. Gugercin, A.C. Antoulas, and C. Beattie
(2008) - SIAM Journal on Matrix Analysis and Applications.
"""
FOM2() = ContinuousSISOTransferFunction(
    [41, 50, 140],
    [1, 11, 111, 110, 100],
)

"""
H2 model reduction for SISO systems,
B. De Moor, P. Van Overschee, and G. Schelfout
(1993) - Proc. of the 12th IFAC World Congress, Sydney, Australia.
"""
FourDisk() = DiscreteSISOTransferFunction(
    [0.0448, 0.2368, 0.0013, 0.0211, 0.2250, 0.0219],
    [1, -1.2024, 2.3675, -2.0039, 2.2337, -1.0420, 0.8513],
)

"""
Globally Optimal H2-Norm Model Reduction: A Numerical Linear Algebra Approach,
O.M.Agudelo, V. Christof, B. De Moor
(2021) - IFAC-PapersOnLine.
"""
Mauricio() = ContinuousSISOTransferFunction(
    [1, 9, -10],
    [1, 12, 49, 78],
)

"""
Rational L_2 Approximation: A Non-Gradient Algorithm,
A. Lepschy, G.A. Mian, G. Pinato, U. Viaro
(1991) - Proc. of the 30th Conference on Decision and Control, Brighton, England.
"""
RationalL2_SISO() = ContinuousSISOTransferFunction(
    [2, 11.5, 57.75, 178.625, 345.5, 323.625, 94.5],
    [1, 10, 46, 130, 239, 280, 194, 60],
)

struct ContinuousStateSpace{T}
    A::Matrix{T}
    B::Matrix{T}
    C::Matrix{T}
    D::Matrix{T}
end

function Base.convert(::Type{ContinuousSISOTransferFunction{T}}, s::ContinuousStateSpace) where {T}
    sys = ControlSystemsBase.tf(ControlSystemsBase.ss(s.A, s.B, s.C, s.D))
    sys = sys.matrix[1]
    return ContinuousSISOTransferFunction{T}(
        reverse(Polynomials.coeffs(sys.num)),
        reverse(Polynomials.coeffs(sys.den)),
    )
end

"""
```julia
# Solution from paper:
A = -4998.078625
B = 100.000194
C = 100.000194
println("Solution proposed in literature:")
println(ControlSystemsBase.tf(ControlSystemsBase.ss(A, B, C, 0)))
```

Homotopy methods for solving the optimal projection equations for the H2 reduced order model problem,
D. Zigic, L.T. Watson, E.G. Collins, D.S. Bernstein
(1992) - International Journal of Control.
"""
HomotopyH2_2() = ContinuousStateSpace(
    [-0.05 -0.99; -0.99 -5000],
    reshape([1.0; 100], 2, 1),
    [1.0 100],
    zeros(1, 1),
)

model_order(tf::SISOTransferFunction) = length(tf.den) - 1

function composeSystemWalsh(tf::SISOTransferFunction, q::Int)
    n = model_order(tf)
    @assert (n >= length(tf.num)) # (strictly) causal model
    @assert q < n

    # Symbolic variables:
    @polyvar a[1:q] b[0:q-1] r[0:n-q-1]
    aa = [1; a[1:end]]

    # Compose coefficient vector of Walsh-polynomial:
    tmp1 = convolveSym(tf.num, aa)
    tmp2 = convolveSym(b, tf.den)
    if tf isa DiscreteSISOTransferFunction
        # Mirror over unit circle (discrete-time)
        tmp3 = convolveSym(convolveSym(reverse(aa), reverse(aa)), r)
    else
        # Flip around imaginary axis (continuous-time)
        aa[end-1:-2:1] = -aa[end-1:-2:1]
        tmp3 = convolveSym(convolveSym(aa, aa), r)
    end
    sys = addPoly(addPoly(tmp1, -tmp2), -tmp3)
    vars = [a; b; r]
    return sys, vars
end


function analyseSolutions(x::Vector{Vector{Float64}}, tf::SISOTransferFunction, q::Int, vars::Vector{PolyVar{true}}, tol::Float64=1e-5)
    n = length(tf.den) - 1
    k = length(x)

    @polyvar s
    b_s, a_s = tf.num' * monomials(s, 0:n-1), tf.den' * monomials(s, 0:n)
    h_exact(y) = b_s(s => y) ./ a_s(s => y)

    fonc = Vector{Vector{ComplexF64}}(undef, k)
    poles = Vector{Vector{ComplexF64}}(undef, k)
    costValues = Vector{Float64}(undef, k)
    validSols = []

    for i = eachindex(x)
        println("----- Solution i: $(i) -----")
        sol = x[i]

        b_hat_s, a_hat_s = sol[q+1:q+q]' * monomials(s, 0:q-1), [1; sol[1:q]]' * monomials(s, 0:q)
        h_approx(y) = b_hat_s(s => y) ./ a_hat_s(s => y)

        # Use Polynomial.jl package to calculate roots of univariate polynomial:
        println("Calculating poles:")
        @time poles[i] = PolynomialRoots.roots(reverse([1; sol[1:q]]))

        # Verify whether FONC are satisfied:
        println("Estimating FONC:")
        tmp1, tmp2 = 0.0, 0.0
        @time begin
            for j = collect(1:q)
                # Optimal zeros:
                if tf isa DiscreteSISOTransferFunction
                    # Use 'atol' because I might be 0.
                    I, est = QuadGK.quadgk(y -> (h_exact(exp(y * im)) - h_approx(exp(y * im))) * conj(1 ./ (exp(y * im) - poles[i][j])), -pi, pi, atol=1e-8)
                else
                    I, est = QuadGK.quadgk(y -> (h_exact(y * im) - h_approx(y * im)) * conj(1 ./ (y * im - poles[i][j])), -Inf, Inf, atol=1e-8)
                end
                tmp1 += abs(1 / (2 * pi) * I)
                # Optimal poles: 
                if tf isa DiscreteSISOTransferFunction
                    I, est = QuadGK.quadgk(y -> (h_exact(exp(y * im)) - h_approx(exp(y * im))) * conj(1 ./ (exp(y * im) - poles[i][j])^2), -pi, pi, atol=1e-8)
                else
                    I, est = QuadGK.quadgk(y -> (h_exact(y * im) - h_approx(y * im)) * conj(1 ./ (y * im - poles[i][j])^2), -Inf, Inf, atol=1e-8)
                end
                tmp2 += abs(1 / (2 * pi) * I)
            end
        end
        fonc[i] = [tmp1; tmp2]
        (sum(abs.(fonc[i])) < tol) ? (push!(validSols, i)) : ()

        println("Computing J:")
        if tf isa DiscreteSISOTransferFunction
            I, est = QuadGK.quadgk(y -> (h_approx(exp(y * im)) - h_exact(exp(y * im))) * conj(h_approx(exp(y * im)) - h_exact(exp(y * im))), -pi, pi, rtol=1e-10)
        else
            I, est = QuadGK.quadgk(y -> (h_approx(y * im) - h_exact(y * im)) * conj(h_approx(y * im) - h_exact(y * im)), -Inf, Inf, rtol=1e-10)
        end
        println(I)
        costValues[i] = sqrt(1 / (2 * pi) * I)
    end

    # Norm of higher order function:
    if tf isa DiscreteSISOTransferFunction
        I, est = QuadGK.quadgk(y -> (h_exact(exp(y * im)) * conj(h_exact(exp(y * im)))), -pi, pi, rtol=1e-8)
    else
        I, est = QuadGK.quadgk(y -> (h_exact(y * im) * conj(h_exact(y * im))), -Inf, Inf, rtol=1e-8)
    end
    h_norm = sqrt(1 / (2 * pi) * I)


    # Print results:
    println("\n")
    if tf isa DiscreteSISOTransferFunction
        println("---- Discrete-time MOR problem ----")
    else
        println("---- Continuous-time MOR problem ----")
    end
    println("** Norm of the exact tf: $(h_norm)")

    res = TypedTables.Table(i=collect(1:k), f=costValues, f_rel=costValues ./ h_norm, opt=fonc, poles=poles, sol=x)

    println("** Minimizing solution: $(findmin(costValues))")

    return res, validSols
end

function convertToSS(sol::Vector{Float64}, q::Int)
    den, num = [1; sol[1:q]], sol[q+1:2q]
    sys = ControlSystemsBase.tf(num, den)
    return ControlSystemsBase.ssdata(ControlSystemsBase.ss(sys))
end


"""
    convolveSym(x::Vector, y::Vector)

Return the coefficients (in decreasing power order) of the
product of polynomials with coefficients `x` and `y`.
(in decreasing power order).
"""
function convolveSym(x::Vector, y::Vector)
    # Calculates the convolution of two vectors (elements can be complex).

    n, m = length(x), length(y)
    # `x` has power `n - 1` and `y` has power `m - 1` so the product
    # has power `n + m - 2`
    l = n + m - 1

    # Pad with zeros:
    # We add zeros for the coefficients of negative powers so that we can use them in convolution
    x = [x; zeros(l - n, 1)]
    y = [y; zeros(l - m, 1)]

    return map(1:l) do i
        # Coefficient of power `s^(l - i)`
        return LinearAlgebra.dot(view(x, 1:i), view(y, i:-1:1))
    end
end

"""
    addPoly(x::Vector, y::Vector)

Return the coefficients (in decreasing power order) of the
product of polynomials with coefficients `x` and `y`.
(in decreasing power order).
"""
function addPoly(x::Vector{Polynomial{true,Float64}}, y::Vector{Polynomial{true,Float64}})
    # Add two polynomials (given as a coefficient vector)

    n, m = length(x), length(y)
    # We add zeros for higher powers
    if m < n
        y = [zeros(n - m, 1); y]
    else
        x = [zeros(m - n, 1); x]
    end
    return vec(x + y)
end


"""
Notes on computation time:

* q = 1-3 feasible (+-20 seconds).
* q = 2 : 500 paths
* q = 3 : 4 000 paths
* q = 4 : 33 000 paths (1 thread: ETA 3h)
* q = 5 : 164 000 paths (6 threads: ETA 4-5h)

Homotopy methods for solving the optimal projection equations for the H2 reduced order model problem,
D. Zigic, L.T. Watson, E.G. Collins, D.S. Bernstein
(1992) - International Journal of Control.
"""
function HomotopyH2_17()
    # System matrix:
    A = zeros(Float64, 17, 17)
    A[1, 1] = -0.031978272
    A[2, 2] = -0.031978272
    A[1, 17] = 0.0097138566
    A[3, 3] = -5.152212
    A[4, 4] = -5.152212
    A[3, 17] = -0.021760771
    A[5, 5] = -0.1351159
    A[6, 6] = -0.1351159
    A[5, 17] = -0.02179972
    A[7, 7] = -0.42811443
    A[8, 8] = -0.42811443
    A[7, 17] = -0.01042631
    A[9, 9] = -0.064896745
    A[10, 10] = -0.064896745
    A[9, 17] = -0.030531575
    A[11, 11] = -0.048520356
    A[12, 12] = -0.048520356
    A[11, 17] = -0.016843335
    A[13, 13] = -0.036781718
    A[14, 14] = -0.036781718
    A[13, 17] = -0.1248007
    A[15, 15] = -0.025112482
    A[16, 16] = -0.025112482
    A[15, 17] = -0.035415526
    A[17, 17] = -92.399784
    A[1, 2] = -78.54
    A[2, 1] = 78.54
    A[2, 17] = -0.0060463517
    A[3, 4] = -51.457677
    A[4, 3] = 51.457677
    A[4, 17] = -0.0054538246
    A[5, 6] = -15.417859
    A[6, 5] = 15.417859
    A[6, 17] = -0.015063913
    A[7, 8] = -14.698408
    A[8, 7] = 14.698408
    A[8, 17] = -0.0088479697
    A[9, 10] = -12.077045
    A[10, 9] = 12.077045
    A[10, 17] = -0.030260987
    A[11, 12] = -8.9654448
    A[12, 11] = 8.9654448
    A[12, 17] = -0.011449591
    A[13, 14] = -4.9057426
    A[14, 13] = 4.9057426
    A[14, 17] = -0.0005136047
    A[15, 16] = -3.8432892
    A[16, 15] = 3.8432892
    A[16, 17] = -0.028115589

    # Input matrix
    B = [1.8631111; -1.1413786; -1.2105758; 0.31424169; 0.013307797; -0.211128913; 0.19552894; -0.037391511; -0.01049736; -0.011486242; -0.029376402; 0.0082391613; -0.012609562; -0.0022040505; -0.030853234; 0.011671662; 0]

    # Output matrix
    C = [-0.0097138566 0.0060463517 0.021760771 -0.0054538246 0.02179972 0.015063913 -0.01042631 -0.0088479697 0.030531575 0.030260987 0.016843335 0.011449591 0.1248007 -0.0005136047 0.035415526 0.028115589 184.79957]
    return ContinuousStateSpace(A, B, C, zeros(1, 1))
end

end
