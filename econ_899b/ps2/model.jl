# Computational Economics
# Professor JF Houde
# Problem set 2
# Alex von Hafften 
# November 22, 2021


# one-dimensional integration
function integrate_1d(f; a = missing, b = missing, method = "quadrature")

    result = 0

    if method == "quadrature"

        # define functions to translate the (0, 1) interval into appropriate interval
        # plus rho_p, the jacobian of the translation 
        rho(u) = u
        rho_p(u) = 1
        
        if !ismissing(a) & !ismissing(b)
            rho(u) = (b - a) * u + a
            rho_p(u) = b - a
        elseif ismissing(a) & ismissing(b)
            rho(u) = -log((1-u)/u)
            rho_prime(u) = -(u/(1-u)) * ((-1)/u - (1-u)/(u^2))
        elseif !ismissing(a) & ismissing(b)
            rho(u) = -log(1 - u) + a
            rho_prime(u) = 1/(1-u)
        elseif ismissing(a) & !ismissing(b)
            rho(u) = log(u) + b
            rho_prime(u) = 1/u
        else
            error("Error in one-dimensional quadrature integration.")
        end

        # sum over grid points
        for i = 1:size(KPU_1d)[1]
            result += KPU_1d[i, :weight] * f(rho(KPU_1d[i, :x1])) * rho_p(KPU_1d[i, :x1])
        end
    end 

    return result

end


# one-dimensional integration
function integrate_2d(f; a_1 = missing, b_1 = missing, a_2 = missing, b_2 = missing, method = "quadrature")

    result = 0

    if method == "quadrature"
        # translation function and jacobian for first dimension
        rho_1(u) = u
        rho_p_1(u) = 1
        if !ismissing(a_1) & !ismissing(b_1)
            rho_1(u) = (b_1 - a_1) * u + a_1
            rho_p_1(u) = b_1 - a_1
        elseif ismissing(a_1) & ismissing(b_1)
            rho_1(u) = -log((1-u)/u)
            rho_p_1(u) = -(u/(1-u)) * ((-1)/u - (1-u)/(u^2))
        elseif !ismissing(a_1) & ismissing(b_1)
            rho_1(u) = -log(1 - u) + a_1
            rho_p_1(u) = 1/(1-u)
        elseif ismissing(a_1) & !ismissing(b_1)
            rho_1(u) = log(u) + b_1
            rho_p_1(u) = 1/u
        else
            error("Error in two-dimensional quadrature integration.")
        end

        # translation function and jacobian for second dimension
        rho_2(u) = u
        rho_p_2(u) = 1
        if !ismissing(a_2) & !ismissing(b_2)
            rho_2(u) = (b_2 - a_2) * u + a_2
            rho_p_2(u) = b_2 - a_2
        elseif ismissing(a_2) & ismissing(b_2)
            rho_2(u) = -log((1-u)/u)
            rho_p_2(u) = -(u/(1-u)) * ((-1)/u - (1-u)/(u^2))
        elseif !ismissing(a_2) & ismissing(b_2)
            rho_2(u) = -log(1 - u) + a_2
            rho_p_2(u) = 1/(1-u)
        elseif ismissing(a_2) & !ismissing(b_2)
            rho_2(u) = log(u) + b_2
            rho_p_2(u) = 1/u
        else
            error("Error in two-dimensional quadrature integration.")
        end

        for i = 1:size(KPU_2d)[1]
            result += KPU_2d[i, :weight] * f(rho_1(KPU_2d[i, :x1]), rho_2(KPU_2d[i, :x2])) * rho_p_1(KPU_2d[i, :x1]) * rho_p_2(KPU_2d[i, :x2])
        end
    end 

    return result

end

# standard normal pdf
function ϕ(x::Float64)
    1/sqrt(2 * π) * exp((-1/2)*x^2)
end

# standard normal cdf
function Φ(x::Float64)
    integrate_1d(ϕ; b = x)
end

# observation-level likelihood
function likelihood(γ::Array{Float64, 1}, β::Array{Float64, 1}, ρ::Float64, α_0::Float64, α_1::Float64, 
    t::Float64, x::Array{Float64, 1}, z::Array{Float64, 1}; method = "quadrature")
    
    result = 0
    σ_0 = 1/(1 - ρ)

    if t == 1
        (-α_0 - x*β' - z*γ')
    elseif t == 2
        result = 2
    elseif t == 3
        result = 3
    elseif t == 4
        result = 4
    else
        error("Invalid value of t.")
    end
    return result
end

# matrix-level likelihood
function likelihood(γ::Array{Float64, 1}, β::Array{Float64, 1}, ρ::Float64, α_0::Float64, α_1::Float64, 
    t::Array{Float64, 1}, x::Array{Float64, 2}, z::Array{Float64, 2}; method = "quadrature")
    N = size(x)[1]

    result = 0

    for i = 1:N
        result += likelihood(γ, β, ρ, α_0, α_1, t[i], x[i,:], z[i,:])
    end

    return result
end