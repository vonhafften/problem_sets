# Computational Economics
# Professor JF Houde
# Problem set 2
# Alex von Hafften 
# November 22, 2021

# compute likelihood for the matrix
function likelihood(γ::Array{Float64, 1}, β::Array{Float64, 1}, ρ::Float64, α_0::Float64, α_1::Float64, α_2::Float64,
    t::Array{Float64, 2}, x::Array{Float64, 2}, z::Array{Float64, 2}; method = "quadrature")
    
    N = size(x)[1]
    
    # uses pmap; on a test with quadrature, this took about 3 minutes for the dataset
    # result = pmap(i-> likelihood(γ, β, ρ, α_0, α_1, α_2, t[i], x[i,:], z[i,:]), i for i = 1:N)
    
    # uses distributed for loop; on a test with quadrature, this took about 2 minutes for the dataset
    result = SharedArray{Float64}(N)

    if method == "quadrature"
        @sync @distributed for i = 1:N
            result[i] = likelihood_quadrature(γ, β, ρ, α_0, α_1, α_2, t[i], x[i,:], z[i,:])
        end

    elseif method == "ghk"
        @sync @distributed for i = 1:N
            result[i] = likelihood_ghk(γ, β, ρ, α_0, α_1, α_2, t[i], x[i,:], z[i,:])
        end

    elseif method == "accept_reject"
        @sync @distributed for i = 1:N
            result[i] = likelihood_accept_reject(γ, β, ρ, α_0, α_1, α_2, t[i], x[i,:], z[i,:])
        end
    end

    return result
end

function log_likelihood(γ::Array{Float64, 1}, β::Array{Float64, 1}, ρ::Float64, α_0::Float64, α_1::Float64, α_2::Float64,
    t::Array{Float64, 2}, x::Array{Float64, 2}, z::Array{Float64, 2}; method = "quadrature")

    sum(log.(likelihood(γ, β, ρ, α_0, α_1, α_2, t, x, z; method = "quadrature")))

end

function log_likelihood(θ::Array{Float64, 1}, t::Array{Float64, 2}, x::Array{Float64, 2}, z::Array{Float64, 2}; method = "quadrature")
    K_x = size(x)[2]
    K_z = size(z)[2]

    γ = θ[1:K_z]
    β = θ[(K_z+1):(K_x + K_z)]
    ρ = θ[K_x + K_z + 1]
    α_0 = θ[K_x + K_z + 2]
    α_1 = θ[K_x + K_z + 3]
    α_2 = θ[K_x + K_z + 4]

    log_likelihood(γ, β, ρ, α_0, α_1, α_2, t, x, z; method = "quadrature")
end 
