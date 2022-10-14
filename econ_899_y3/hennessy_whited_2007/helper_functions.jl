# Code to solve Hennessy and Whited 2007

# Alex von Hafften
# Sept. 27, 2022

using Parameters

include("structures.jl")

# corporate tax bill
function T_C(x::Float64, P::Primitives)
    if x <= 0.0
        return P.τ_c_n  * x
    else
        return P.τ_c_p * x
    end
end

# cash distribution tax bill
function T_d(x::Float64, P::Primitives)
    if x <= 0.0
        return 0.0
    else
        return P.τ_d_bar/P.ϕ * (x*P.ϕ + exp(-P.ϕ*x) - 1.0)
    end
end

# cost of external equity
function Λ(x::Float64, P::Primitives)
    if x <= 0.0
        return 0.0
    else
        return P.λ_0 + P.λ_1 * x + P.λ_2 * x^2
    end
end
