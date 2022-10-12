# Code to solve Hennessy and Whited 2007

# Alex von Hafften
# Sept. 27, 2022

using Parameters

include("structures.jl")

# corporate tax bill
function T_C(x::Float64, P::Primitives)
    if x <= 0
        return P.τ_c_n  * x
    else
        return P.τ_c_p * x
    end
end

# cash distribution tax bill
function T_d(x::Float64, P::Primitives)
    if x <= 0
        return 0
    else
        return P.τ_d_bar/P.ϕ * (x*P.ϕ + exp(-P.ϕ*x) - 1)
    end
end

# cost of external equity
function Λ(x::Float64, P::Primitives)
    if x <= 0
        return 0
    else
        return P.λ_0 + P.λ_1 * x + P.λ_2 * x^2
    end
end

# compute recovery in default
function compute_R(k_p::Float64, z_p::Float64, w_bar::Float64, P::Primitives)
    temp = z_p * compute_π(k_p, P) - P.δ * k_p
    (1 - P.ξ) * (1 - P.δ) * k_p - z_p * compute_π(k_p, P) - (P.τ_c_p*(temp > 0) + P.τ_c_n*(temp <= 0)) * temp - w_bar
end