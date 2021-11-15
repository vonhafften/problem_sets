# Computational Economics
# Professor JF Houde
# Problem set 2
# Alex von Hafften 
# November 22, 2021

# standard normal pdf
@everywhere function ϕ(x::Float64)
    1/sqrt(2 * π) * exp((-1/2)*x^2)
end

# standard normal cdf
@everywhere function Φ(x::Float64)
    integrate_1d(ϕ, x)
end

# observation-level likelihood
@everywhere function likelihood(γ::Array{Float64, 1}, β::Array{Float64, 1}, ρ::Float64, α_0::Float64, α_1::Float64, α_2::Float64, 
    t::Float64, x::Array{Float64, 1}, z::Array{Float64, 1})
    
    result = 0.0
    σ_0 = 1/(1 - ρ)

    f_2(ε_0)      = Φ(-α_1 - x'*β - z'*γ - ρ*ε_0) * ϕ(ε_0/σ_0) / σ_0
    f_3(ε_0, ε_1) = Φ(-α_2 - x'*β - z'*γ - ρ*ε_1) * ϕ(ε_1 - ρ*ε_0) * ϕ(ε_0/σ_0) / σ_0
    f_4(ε_0, ε_1) = Φ( α_2 + x'*β + z'*γ - ρ*ε_1) * ϕ(ε_1 - ρ*ε_0) * ϕ(ε_0/σ_0) / σ_0

    if t == 1.0
        result = Φ((-α_0 - x'*β - z'*γ)/σ_0)
    elseif t == 2.0
        result = integrate_1d(f_2, α_0 + x'*β + z'*γ)
    elseif t == 3.0
        result = integrate_2d(f_3, α_0 + x'*β + z'*γ, α_1 + x'*β + z'*γ)
    elseif t == 4.0
        # getting probabilities larger than one.
        result = integrate_2d(f_4, α_0 + x'*β + z'*γ, α_1 + x'*β + z'*γ)

        # result_1 = Φ((-α_0 - x'*β - z'*γ)/σ_0)
        # f_2(ε_0) = Φ(-α_1 - x'*β - z'*γ - ρ*ε_0) * ϕ(ε_0/σ_0) / σ_0
        # result_2 = integrate_1d(f_2, α_0 + x'*β + z'*γ)
        # f_3(ε_0, ε_1) = Φ(-α_2 - x'*β - z'*γ - ρ*ε_1) * ϕ(ε_1 - ρ*ε_0) * ϕ(ε_0/σ_0) / σ_0
        # result_3 = integrate_2d(f_3, α_0 + x'*β + z'*γ, α_1 + x'*β + z'*γ)
        # result = 1 - result_1 - result_2 - result_3
    else
        error("Invalid value of t.")
    end
    return result
end

# compute likelihood for the matrix
function likelihood(γ::Array{Float64, 1}, β::Array{Float64, 1}, ρ::Float64, α_0::Float64, α_1::Float64, α_2::Float64,
    t::Array{Float64, 2}, x::Array{Float64, 2}, z::Array{Float64, 2})
    
    N = size(x)[1]
    
    # uses pmap; on a test with quadrature, this took about 3 minutes for the dataset
    # result = pmap(i-> likelihood(γ, β, ρ, α_0, α_1, α_2, t[i], x[i,:], z[i,:]), i for i = 1:N)
    
    # uses distributed for loop; on a test with quadrature, this took about 2 minutes for the dataset
    result = SharedArray{Float64}(N)
    @sync @distributed for i = 1:N
        result[i] = likelihood(γ, β, ρ, α_0, α_1, α_2, t[i], x[i,:], z[i,:])
    end

    return result
end

function log_likelihood(γ::Array{Float64, 1}, β::Array{Float64, 1}, ρ::Float64, α_0::Float64, α_1::Float64, α_2::Float64,
    t::Array{Float64, 2}, x::Array{Float64, 2}, z::Array{Float64, 2})

    sum(log.(likelihood(γ, β, ρ, α_0, α_1, α_2, t, x, z)))

end

function log_likelihood(θ::Array{Float64, 1}, t::Array{Float64, 2}, x::Array{Float64, 2}, z::Array{Float64, 2})
    K_x = size(x)[2]
    K_z = size(z)[2]

    γ = θ[1:K_z]
    β = θ[(K_z+1):(K_x + K_z)]
    ρ = θ[K_x + K_z + 1]
    α_0 = θ[K_x + K_z + 2]
    α_1 = θ[K_x + K_z + 3]
    α_2 = θ[K_x + K_z + 4]

    log_likelihood(γ, β, ρ, α_0, α_1, α_2, t, x, z)
end 
