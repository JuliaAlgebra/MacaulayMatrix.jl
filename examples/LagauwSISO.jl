# ----------------------------- SISO --------------------------------
# Let's try an example from the SISO misfit modeling problem.
# The solution set of this system of equations is positive dimensional over C. 

using LinearAlgebra
using TypedPolynomials
using MacaulayMatrix
using JuMP
using MultivariateMoments
using DynamicPolynomials


function composeSystemRiSVD(y::Vector{Float64}, u::Vector{Float64}, n::Integer)
    @assert length(y) == length(u)
    N = length(y)

    DynamicPolynomials.@polyvar a[2:n+1] b[1:n+1]
    a = reverse(a)
    a = [a; 1]
    b = reverse(b)

    Ta_Nn = zeros(eltype(a), (N - n), N)
    Tb_Nn = zeros(eltype(a), (N - n), N)
    for i = axes(Ta_Nn, 1), j in eachindex(a)
        if i + j - 1 in axes(Ta_Nn, 2)
            Ta_Nn[i, i+j-1] = a[j]
            Tb_Nn[i, i+j-1] = b[j]
        end
    end

    DynamicPolynomials.@polyvar l[1:N-n]
    Tl_n1 = zeros(eltype(a), n + 1, N)
    for i = axes(Tl_n1, 1), j in eachindex(l)
        if i + j - 1 in axes(Tl_n1, 2)
            Tl_n1[i, i+j-1] = l[j]
        end
    end

    fonc = Vector{polynomial_type(eltype(a), Float64)}(undef, N + n + 3)
    fonc[1:N-n] = (Ta_Nn * transpose(Ta_Nn) + Tb_Nn * transpose(Tb_Nn)) * l - Ta_Nn * y + Tb_Nn * u
    fonc[N-n+1:N+1] = Tl_n1 * y - Tl_n1 * transpose(Tl_n1) * a
    fonc[N+2:N+n+2] = Tl_n1 * u + Tl_n1 * transpose(Tl_n1) * b
    fonc[N+n+3] = transpose(a) * a + transpose(b) * b - 1 # Normalization

    return fonc
end

y = [1, 2, 3, 4, 5]
u = [5, -4, 3, -2, 1]
y = convert(Vector{Float64}, y)
u = convert(Vector{Float64}, u)
n = 1

sysRiSVD = composeSystemRiSVD(y, u, n)

# Try to solve using Macaulay:
sols = solve_system(sysRiSVD, column_maxdegree=10, print_level=3)