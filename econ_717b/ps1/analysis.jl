# Alex von Hafften
# Problem set 1
# ECON 717B: Applied Econometrics
# Professor Matt Wiswall
# April 7, 2022

using Parameters, Distributions

@with_kw struct parameters_0
    π_1 = 1.0
    π_2 = 1.0
    μ_1 = 1.0
    μ_2 = 1.0
    σ_1 = 1.0
    σ_2 = 1.0
    ρ   = 0.5
end

mutable struct data
    N::Int64             # number of trials
    ε_1::Vector{Float64} # idiosyncratic component for skill in occupation 1
    ε_2::Vector{Float64} # idiosyncratic component for skill in occupation 2
    S_1::Vector{Float64} # skill in occupation 1
    S_2::Vector{Float64} # skill in occupation 2
    W_1::Vector{Float64} # potential wages in occupation 1
    W_2::Vector{Float64} # potential wages in occupation 2
    D::Vector{Float64}   # occupation choice (observed)
    W::Vector{Float64}   # wages (observed)
end

function simulate_date()
    p = parameters_0()

    σ_11 = σ_1 * σ_1
    σ_22 = σ_2 * σ_2
    σ_12 = ρ * σ_1 * σ_2

    N = 1000

    ε_distribution = MvNormal([p.μ_1, p.μ_2], [σ_11 σ_12, σ_12 σ_22])

end