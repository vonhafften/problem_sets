# Computational Economics
# Professor JF Houde
# Problem set 2
# Alex von Hafften 
# November 22, 2021

# generates halton sequence of length n using base (base must be prime)
function halton(base::Int64, n::Int64)
    
    m, d = 0, 1
    burn_in = 200
    result = zeros(burn_in + n)

    # Adapted from code here: https://en.wikipedia.org/wiki/Halton_sequence
    for i = 1:(burn_in + n)
        x = d - m
        if x == 1
            m = 1
            d *= base
        else
            y = d / base
            while x <= y
                y /= base
            end
            m = (base + 1) * y - x
        end
        result[i] = m / d
    end

    result[(burn_in + 1):(burn_in + n)]

end

# returns three iid uniform(0, 1) RVs using halton for GHk
function initialize_ghk(;use_halton = true)
    n_trials = 100

    # pulls independent uniform shocks
    if use_halton
        
        u_0 = halton(5,  n_trials)
        u_1 = halton(7,  n_trials)
        u_2 = halton(11, n_trials)

    else

        uniform_distibution = Uniform(0, 1)
        
        Random.seed!(1)
        u_0 = rand(uniform_distibution, n_trials)
        
        Random.seed!(2)
        u_1 = rand(uniform_distibution, n_trials)

        Random.seed!(3)
        u_2 = rand(uniform_distibution, n_trials)
    
    end

    return [u_0, u_1, u_2]
end

# returns correlated ε for accept-reject
function initialize_accept_reject(ρ::Float64; use_halton = true)

    u_0, u_1, u_2 = initialize_ghk(;use_halton)

    n_trials = length(u_0)

    # initialize vectors to store normal shocks
    η_0 = zeros(n_trials)
    η_1 = zeros(n_trials)
    η_2 = zeros(n_trials)

    # uses Φ_inverse function to transform from uniform to normal
    for i in 1:n_trials
        η_0[i] = Φ_inverse(u_0[i])
        η_1[i] = Φ_inverse(u_1[i])
        η_2[i] = Φ_inverse(u_2[i])
    end

    # Define correlated errors
    ε_0 = η_0 .* (1/(1 - ρ)^2)
    ε_1 = ρ .* ε_0 .+ η_1
    ε_2 = ρ .* ε_1 .+ η_2

    return [ε_0, ε_1, ε_2]
end

# reads in grid points for quadrature integration
function initialize_quadrature_integration()

    # quadrature nodes and weights
    KPU_1d = DataFrame(CSV.File("PS2/KPU_d1_l20.csv"))
    KPU_2d = DataFrame(CSV.File("PS2/KPU_d2_l20.csv"))

    return [KPU_1d, KPU_2d]
end

# compute likelihood for the matrix
function likelihood(α_0::Float64, α_1::Float64, α_2::Float64,  β::Array{Float64, 1}, γ::Float64, ρ::Float64, 
    t::Array{Float64, 2}, x::Array{Float64, 2}, z::Array{Float64, 2},  KPU_1d, KPU_2d,
    u_0::Array{Float64, 1}, u_1::Array{Float64, 1}, u_2::Array{Float64, 1},
    ε_0::Array{Float64, 1}, ε_1::Array{Float64, 1}, ε_2::Array{Float64, 1}; method = "quadrature")
    
    N = size(x)[1]
    
    # uses distributed for loop; on a test with quadrature, this took about 2 minutes for the dataset
    result = SharedArray{Float64}(N)

    if method == "quadrature"
        println("Evaluating likelihoods using quadrature integration method...")

        @showprogress @distributed for i = 1:N
            result[i] = likelihood_quadrature(α_0, α_1, α_2, β, γ, ρ, t[i], x[i,:], z[i,:], KPU_1d, KPU_2d)
        end

    elseif method == "ghk"
        println("Evaluating likelihoods using GHK method...")

        @showprogress @distributed for i = 1:N
            result[i] = likelihood_ghk(α_0, α_1, α_2, β, γ, ρ, t[i], x[i,:], z[i,:], u_0, u_1, u_2)
        end

    elseif method == "accept_reject"
        println("Evaluating likelihoods using accept-reject method...")

        @showprogress @distributed for i = 1:N
            result[i] = likelihood_accept_reject(α_0, α_1, α_2, β, γ, ρ, t[i], x[i,:], z[i,:], ε_0, ε_1, ε_2)
        end
    else
        error("Specify valid method.")
    end

    return result
end

function log_likelihood(θ::Array{Float64, 1}, t::Array{Float64, 2}, x::Array{Float64, 2}, z::Array{Float64, 2},  KPU_1d, KPU_2d,
    u_0::Array{Float64, 1}, u_1::Array{Float64, 1}, u_2::Array{Float64, 1},
    ε_0::Array{Float64, 1}, ε_1::Array{Float64, 1}, ε_2::Array{Float64, 1}; method = "quadrature")

    K_x = size(x)[2]

    α_0 = θ[1]
    α_1 = θ[2]
    α_2 = θ[3]
    β   = θ[4:(K_x+3)]
    γ   = θ[K_x+4]
    ρ   = θ[K_x+5]

    ll = sum(log.(likelihood(α_0, α_1, α_2, β, γ, ρ, t, x, z, KPU_1d, KPU_2d, u_0, u_1, u_2, ε_0, ε_1, ε_2; method = method)))

    println("Log-likelihood = ", ll)

    println("α_0 = ", α_0)
    println("α_1 = ", α_1)
    println("α_2 = ", α_2)
    println("β = ", β)
    println("γ = ", γ)
    println("ρ = ", ρ)

    println("******************************************************")

    return ll
end 
