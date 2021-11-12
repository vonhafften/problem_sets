# Computational Economics
# Professor JF Houde
# Problem set 2
# Alex von Hafften 
# November 22, 2021

# standard normal pdf
function ϕ(x::Float64)
    1/sqrt(2 * π) * exp((-1/2)*x^2)
end

# standard normal cdf
function Φ(x::Float64)
    integrate_1d(ϕ; b = x)
end

# observation-level likelihood
function likelihood(γ::Array{Float64, 1}, β::Array{Float64, 1}, ρ::Float64, α_0::Float64, α_1::Float64, α_2::Float64, 
    t::Float64, x::Array{Float64, 1}, z::Array{Float64, 1}; method = "quadrature")
    
    result = 0.0
    σ_0 = 1/(1 - ρ)

    if t == 1.0
        result = Φ((-α_0 - x'*β - z'*γ)/σ_0)
    elseif t == 2.0
        f_2(ε_0) = Φ(-α_1 - x'*β - z'*γ - ρ*ε_0) * ϕ(ε_0/σ_0) / σ_0
        result = integrate_1d(f_2, b = α_0 + x'*β + z'*γ)
    elseif t == 3.0
        f_3(ε_0, ε_1) = Φ(-α_2 - x'*β - z'*γ - ρ*ε_1) * ϕ(ε_1 - ρ*ε_0) * ϕ(ε_0/σ_0) / σ_0
        result = integrate_2d(f_3; b_0 = α_0 + x'*β + z'*γ, b_1 = α_1 + x'*β + z'*γ)
    elseif t == 4.0
        f_4(ε_0, ε_1) = Φ(α_2 + x'*β + z'*γ - ρ*ε_1)  * ϕ(ε_1 - ρ*ε_0) * ϕ(ε_0/σ_0) / σ_0
        result = integrate_2d(f_4; b_0 = α_0 + x'*β + z'*γ, b_1 = α_1 + x'*β + z'*γ)
    else
        error("Invalid value of t.")
    end
    return result
end

# compute likelihood
function log_likelihood(γ::Array{Float64, 1}, β::Array{Float64, 1}, ρ::Float64, α_0::Float64, α_1::Float64, α_2::Float64,
    t::Array{Float64, 2}, x::Array{Float64, 2}, z::Array{Float64, 2}; method = "quadrature")
    N = size(x)[1]

    result = 0

    for i = 1:N
        # println(i)
        result += log(likelihood(γ, β, ρ, α_0, α_1, α_2, t[i], x[i,:], z[i,:]))
    end

    return result
end

function log_likelihood(θ::Array{Float64, 1}, t::Array{Float64, 2}, x::Array{Float64, 2}, z::Array{Float64, 2})

end