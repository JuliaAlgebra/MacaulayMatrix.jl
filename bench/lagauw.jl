using LinearAlgebra
using DynamicPolynomials

struct Lagauw1
    mitigate_rational::Bool
end

function system(bench::Lagauw1)
    # Given data (w = [u, y]):
    w = [
        1  2 3  4 5
        5 -4 3 -2 1
    ]
    N = size(w, 2)

    # Reshape to w = [y1, u1, y2, u2, ...]
    w = [w[2, :]; w[1, :]]

    # Define variables:
    # @PolyVar (DynamicPolynomials) does not allow to compute with symbolic matrices
    # @var (HomotopyContinuation)
    # @variables (Symbolics) 
    n = 1
    @polyvar a[1:n] b[0:n] # a[1] = a_1, b[1] = b_0, b[2] = b_1
    vars = [a; b[2]; 1; b[1]] # [a_1, b_1, 1, b_0]

    # Create system:
    # 1. Toeplitz matrix
    TermType = termtype(a[1], Float64)
    PolynomialType = polynomialtype(a[1], Float64)
    T = zeros(TermType, 2 * (N - n), 2 * N)
    for i = axes(T, 1), j in eachindex(vars)
        #if i + j - 1 <= 2 * N
        if i + j - 1 in axes(T, 2)
            T[i, i+j-1] = vars[j]
        end
    end
    T = T[1:2:size(T, 1), :] # Select 'double-shifts' only...
    @assert (size(T, 1) == N - n && size(T, 2) == 2 * N) "Size of T incorrect..."

    D = T * transpose(T)

    if bench.mitigate_rational
        @polyvar f[1:N-n] fii[1:(2*n+1)*(N-n)]
        fi = zeros(TermType, 2 * n + 1, N - n)
        for i = 1:2*n+1
            fi[i, :] = fii[1+(i-1)*(N-n):i*(N-n)]
        end
        # Each row of fi contains the part. derivs. of f wrt. a model variable (a_i or b_i)
        # fi = [fa1^T;
        #       fa2^T;
        #       fb0^T;
        #       fb1^T;
        #       fb2^T]

        # --- Building the system ---
        system = Vector{PolynomialType}(undef, (2 * n + 1) * (N - n + 1) + (N - n))

        # 1) N-n eqs to fix 'f':
        # (defining the auxiliary variable f)
        system[1:N-n] = D * f - T * w

        # 2) First line of TK85 eq. (11): (N-n) eqs per model variable
        # (defining the partial derivs of f wrt. the model variables)
        # a_i
        for i = 1:n
            system[N-n+1+(i-1)*(N-n):N-n+i*(N-n)] = differentiate.(D, a[i]) * f - differentiate.(T, a[i]) * w + D * fi[i, :]
        end
        # b_i
        start = N - n + n * (N - n)
        for i = 1:n+1
            system[start+1+(i-1)*(N-n):start+i*(N-n)] = differentiate.(D, b[i]) * f - differentiate.(T, b[i]) * w + D * fi[n+i, :]
        end

        # 3) Second line of TK85 eq. (11): 1 eq per model variable
        # (The FONC as a function of model variables and f and fi)
        start = N - n + n * (N - n) + (n + 1) * (N - n)
        for i = 1:n
            tmp = transpose(w) * transpose(differentiate.(T, a[i])) * f + transpose(w) * transpose(T) * fi[i, :]
            system[start+i] = tmp[1]
        end
        for i = 1:n+1
            tmp = transpose(w) * transpose(differentiate.(T, b[i])) * f + transpose(w) * transpose(T) * fi[n+i, :]
            system[start+n+i] = tmp[1]
        end
    else
        # Tried this and it failed:
        # 10seconds to compute admissible tuple
        # constructing Sylvester matrix: 7.627104 seconds (70.58 M allocations: 7.560 GiB, 12.40% gc time, 0.10% compilation time), Sylv is a matrix of size (18238, 17730)
        # computing cokernel: 1836.309116 seconds (13 allocations: 14.258 GiB, 0.04% gc time), rank = 12830, relative size of last nz singular value = 3.48711582338338e-10, gap = 1.0581115621040016e6, N has size (5408, 18238)
        # Checking criterion: !!!!!!!!! criterion violated !!!!!!!!! Criterion was not satisfied. 2037.464780 seconds (100.48 M allocations: 37.539 GiB, 0.12% gc time, 0.38% compilation time)
        # (nothing, [0 0 0; 0 0 1; 0 1 0; 1 0 0], [[0 1 0; 0 1 1; … ; 43 0 1; 43 0 2], [0 1 0; 0 1 1; … ; 30 0 0; 30 0 1], [0 1 0; 0 1 1; … ; 30 0 0; 30 0 1], [0 0 0; 0 0 1; … ; 29 0 1; 29 0 2]], [0 1 0; 0 1 1; … ; 44 0 1; 44 0 2])


        # 2. Cofactor matrix:
        C = zeros(PolynomialType, N - n, N - n)
        for i in 1:N-n, j in 1:N-n
            C[i, j] = (-1)^(i + j) * det(D[union(1:i-1, i+1:N-n), union(1:j-1, j+1:N-n)])
        end
        # Adjugate matrix is transpose of cofactor matrix
        C = transpose(C)

        # 3. FONC:
        # We assume det(D) \neq 0
        system = map(1:(2n + 1)) do i
            v = i <= n ? a[i] : b[i - n]
            2 * transpose(w) * transpose(T) * C * differentiate.(T, v) * w * det(D) - transpose(w) * transpose(T) * C * differentiate.(D, v) * C * T * w
        end

        # Require det(D)*f = 1 (with an auxiliary variable f), making sure that the det is not equal to zero, hopefully reducing the positive dimensionality.
        @polyvar f
        push!(system, det(D) * f - 1)
    end
    return system
end

s0 = system(Lagauw1(false))
s1 = system(Lagauw1(true))

using Macaulay
using MultivariateMoments
solve_system(s1, 4)

import SCS
ν13scs = moment_matrix(s1, SCS.Optimizer, 3)
ν14scs = moment_matrix(s1, SCS.Optimizer, 4)

using MosekTools
ν13mosek = moment_matrix(s1, Mosek.Optimizer, 3)
ν14mosek = moment_matrix(s1, Mosek.Optimizer, 4)
