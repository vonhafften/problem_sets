# Code to solve Hennessy and Whited 2007

# Alex von Hafften
# Sept. 27, 2022

using Parameters

include("structures.jl")

# operating profits
function compute_π(k::Float64, P::Primitives)
    k^P.α
end

# inverse of operating profits derivative
function compute_π_d_inv(profits::Float64, P::Primitives)
    (profits/P.α)^(1/(P.α-1))
end

# corporate tax bill
function compute_T_C(k::Float64, b::Float64, z::Float64, r_tilde::Float64, P::Primitives)
    temp = z * compute_π(k, P) - P.δ * k - r_tilde * b
    if temp <= 0
        return P.τ_c_n  * temp
    else
        return P.τ_c_p * temp
    end
end

# cash distribution tax bill
function compute_T_d(x::Float64, P::Primitives)
    if x <= 0
        return 0
    else
        return P.τ_d_bar * (x + exp(-P.ϕ*x)/P.ϕ - 1/P.ϕ)
    end
end

# net worth
function compute_nw(k_p::Float64, b_p::Float64, z::Float64, z_p::Float64, r_tilde::Float64, P::Primitives)
    return z_p * compute_π(k_p, P) + (1-P.δ)*k_p - compute_T_C(k_p, b_p, z, r_tilde, P) - (1+r_tilde)*b_p
end

# cost of external equity
function compute_Λ(x::Float64, P::Primitives)
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