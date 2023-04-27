using DynamicPolynomials
using LinearAlgebra
import QuadGK
import PolynomialRoots
import TypedTables
import ControlSystemsBase
import Polynomials

function selectProblem(problem::String)
    # Paper: 
    if problem == "F1"
        # Paper: Second order H2 optimal approximation of linear dynamical systems, M.I. Ahmad, M. Frangos, I.M. Jaimoukha (2011) - Proc. of the 18th IFAC World Congress, Milano, Italy.
        b = [-2.9239, -39.5525, -97.5270, -147.1508]
        a = [1, 11.9584, 43.9119, 73.6759, 44.3821]
        isDiscrete = false
    elseif problem == "F2"
        # Paper: Second order H2 optimal approximation of linear dynamical systems, M.I. Ahmad, M. Frangos, I.M. Jaimoukha (2011) - Proc. of the 18th IFAC World Congress, Milano, Italy.
        b = [-1.2805, -6.2266, -12.8095, -9.3373]
        a = [1, 3.1855, 8.9263, 12.2936, 3.1987]
        isDiscrete = false
    elseif problem == "F3"
        # Paper: Second order H2 optimal approximation of linear dynamical systems, M.I. Ahmad, M. Frangos, I.M. Jaimoukha (2011) - Proc. of the 18th IFAC World Congress, Milano, Italy.
        b = [-1.3369, -4.8341, -47.5819, -42.7285]
        a = [1, 17.0728, 84.9908, 122.4400, 59.9309]
        isDiscrete = false
    elseif problem == "Spanos-1992-1"
        # Paper: A New Algorithm for L_2 Optimal Model Reduction, J.T. Spanos, M.H. Milman, D.L. Mingori (1992) - Automatica.
        b = [1, 15, 50]
        a = [1, 5, 33, 79, 50]
        isDiscrete = false
    elseif problem == "F6"
        # Paper: Second order H2 optimal approximation of linear dynamical systems, M.I. Ahmad, M. Frangos, I.M. Jaimoukha (2011) - Proc. of the 18th IFAC World Congress, Milano, Italy.
        b = [41, 50, 140]
        a = [1, 11, 111, 110, 100]
        isDiscrete = false
    elseif problem == "FOM2"
        # Paper: H2 Model Reduction for Large-Scale Linear Dynamical Systems, S. Gugercin, A.C. Antoulas, and C. Beattie (2008) - SIAM Journal on Matrix Analysis and Applications.
        b = [2, 11.5, 57.75, 178.625, 345.5, 323.625, 94.5]
        a = [1, 10, 46, 130, 239, 280, 194, 60]
        isDiscrete = false
    elseif problem == "Spanos-1992-3"
        # Paper: A New Algorithm for L_2 Optimal Model Reduction, J.T. Spanos, M.H. Milman, D.L. Mingori (1992) - Automatica.
        b = [0, 0, 0, 0.00001, 0.0110, 1]
        a = [1, 0.2220, 22.1242, 3.5445, 122.4433, 11.3231, 11.1100]
        isDiscrete = false
    elseif problem == "FourDisk"
        # Paper: H2 model reduction for SISO systems, B. De Moor, P. Van Overschee, and G. Schelfout (1993) - Proc. of the 12th IFAC World Congress, Sydney, Australia.
        b = [0.0448, 0.2368, 0.0013, 0.0211, 0.2250, 0.0219]
        a = [1, -1.2024, 2.3675, -2.0039, 2.2337, -1.0420, 0.8513]
        isDiscrete = true
    elseif problem == "Mauricio"
        # Paper: Globally Optimal H2-Norm Model Reduction: A Numerical Linear Algebra Approach, O.M.Agudelo, V. Christof, B. De Moor (2021) - IFAC-PapersOnLine.
        b = [1, 9, -10]
        a = [1, 12, 49, 78]
        isDiscrete = false
    elseif problem == "HomotopyH2_2"
        # Paper: Homotopy methods for solving the optimal projection equations for the H2 reduced order model problem, D. Zigic, L.T. Watson, E.G. Collins, D.S. Bernstein (1992) - International Journal of Control.
        isDiscrete = false
        A = [-0.05 -0.99; -0.99 -5000]
        B = [1; 100]
        C = [1 100]
        D = 0
        sys = ControlSystemsBase.tf(ControlSystemsBase.ss(A, B, C, D))
        sys = sys.matrix[1]
        b = reverse(Polynomials.coeffs(sys.num))
        a = reverse(Polynomials.coeffs(sys.den))

        # Solution from paper:
        A = -4998.078625
        B = 100.000194
        C = 100.000194
        println("Solution proposed in literature:")
        println(ControlSystemsBase.tf(ControlSystemsBase.ss(A, B, C, 0)))
    elseif problem == "HomotopyH2_17"
        # Paper: Homotopy methods for solving the optimal projection equations for the H2 reduced order model problem, D. Zigic, L.T. Watson, E.G. Collins, D.S. Bernstein (1992) - International Journal of Control.
        isDiscrete = false
        A, B, C = getLargeSystem("HomotopyH2_17")
        sys = ControlSystemsBase.tf(ControlSystemsBase.ss(A, B, C, 0))
        sys = sys.matrix[1]
        b = reverse(Polynomials.coeffs(sys.num))
        a = reverse(Polynomials.coeffs(sys.den))

        # Notes on computation time:
        # q = 1-3 feasible (+-20 seconds).
        # q = 2 : 500 paths
        # q = 3 : 4 000 paths
        # q = 4 : 33 000 paths (1 thread: ETA 3h)
        # q = 5 : 164 000 paths (6 threads: ETA 4-5h)

    elseif problem == "RationalL2_SISO"
        # Paper: Rational L_2 Approximation: A Non-Gradient Algorithm, A. Lepschy, G.A. Mian, G. Pinato, U. Viaro (1991) - Proc. of the 30th Conference on Decision and Control, Brighton, England.
        isDiscrete = false
        b = [2, 11.5, 57.75, 178.625, 345.5, 323.625, 94.5]
        a = [1, 10, 46, 130, 239, 280, 194, 60]
    end

    println("Selected problem: $(problem) (discrete: $(isDiscrete))")
    println("\n")

    return Dict("a" => convert(Vector{Polynomial{true,Float64}}, a), "b" => convert(Vector{Polynomial{true,Float64}}, b), "isDiscrete" => isDiscrete)
end

function composeSystemWalsh(tf::Dict, q::Int)

    n = length(tf["a"]) - 1 # Model order
    @assert (n >= length(tf["b"])) # (strictly) causal model
    @assert q < n

    # Symbolic variables:
    @polyvar a[1:q] b[0:q-1] r[0:n-q-1]
    aa = convert(Vector{Polynomial{true,Float64}}, [1; a[1:end]])
    bb = convert(Vector{Polynomial{true,Float64}}, b)
    rr = convert(Vector{Polynomial{true,Float64}}, r)

    # Compose coefficient vector of Walsh-polynomial:
    tmp1 = convolveSym(tf["b"], aa)
    tmp2 = convolveSym(bb, tf["a"])
    if tf["isDiscrete"]
        # Mirror over unit circle (discrete-time)
        tmp3 = convolveSym(convolveSym(reverse(aa), reverse(aa)), rr)
    else
        # Flip around imaginary axis (continuous-time)
        aa[end-1:-2:1] = -aa[end-1:-2:1]
        tmp3 = convolveSym(convolveSym(aa, aa), rr)
    end
    sys = addPoly(addPoly(tmp1, -tmp2), -tmp3)
    vars = [a; b; r]
    return sys, vars
end


function analyseSolutions(x::Vector{Vector{Float64}}, tf::Dict, q::Int, vars::Vector{PolyVar{true}}, tol::Float64=1e-5)

    n = length(tf["a"]) - 1
    k = length(x)

    @polyvar s
    b_s, a_s = tf["b"]' * monomials(s, 0:n-1), tf["a"]' * monomials(s, 0:n)
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
                if tf["isDiscrete"]
                    # Use 'atol' because I might be 0.
                    I, est = QuadGK.quadgk(y -> (h_exact(exp(y * im)) - h_approx(exp(y * im))) * conj(1 ./ (exp(y * im) - poles[i][j])), -pi, pi, atol=1e-8)
                else
                    I, est = QuadGK.quadgk(y -> (h_exact(y * im) - h_approx(y * im)) * conj(1 ./ (y * im - poles[i][j])), -Inf, Inf, atol=1e-8)
                end
                tmp1 += abs(1 / (2 * pi) * I)
                # Optimal poles: 
                if tf["isDiscrete"]
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
        if tf["isDiscrete"]
            I, est = QuadGK.quadgk(y -> (h_approx(exp(y * im)) - h_exact(exp(y * im))) * conj(h_approx(exp(y * im)) - h_exact(exp(y * im))), -pi, pi, rtol=1e-10)
        else
            I, est = QuadGK.quadgk(y -> (h_approx(y * im) - h_exact(y * im)) * conj(h_approx(y * im) - h_exact(y * im)), -Inf, Inf, rtol=1e-10)
        end
        println(I)
        costValues[i] = sqrt(1 / (2 * pi) * I)
    end

    # Norm of higher order function:
    if tf["isDiscrete"]
        I, est = QuadGK.quadgk(y -> (h_exact(exp(y * im)) * conj(h_exact(exp(y * im)))), -pi, pi, rtol=1e-8)
    else
        I, est = QuadGK.quadgk(y -> (h_exact(y * im) * conj(h_exact(y * im))), -Inf, Inf, rtol=1e-8)
    end
    h_norm = sqrt(1 / (2 * pi) * I)


    # Print results:
    println("\n")
    if tf["isDiscrete"]
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


function convolveSym(x::Vector{Polynomial{true,Float64}}, y::Vector{Polynomial{true,Float64}})
    # Calculates the convolution of two vectors (elements can be complex).

    n, m = length(x), length(y)
    l = n + m - 1

    # Pad with zeros:
    x = [x; zeros(l - n, 1)]
    y = [y; zeros(l - m, 1)]

    res = Array{Polynomial{true,Float64},1}(undef, l)
    for i = 1:l
        tmp = 0
        for j = 1:i
            tmp = tmp + x[j] * y[i+1-j]
        end
        res[i] = tmp
    end
    return res
end

function addPoly(x::Vector{Polynomial{true,Float64}}, y::Vector{Polynomial{true,Float64}})
    # Add two polynomials (given as a coefficient vector)

    n, m = length(x), length(y)
    if m < n
        y = [zeros(Polynomial{true,Float64}, n - m, 1); y]
    else
        x = [zeros(Polynomial{true,Float64}, m - n, 1); x]
    end
    res = vec(x + y)
end



function getLargeSystem(system::String)
    if system == "HomotopyH2_17"
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
    end
    return A, B, C
end
