# Problem Set 3
# Computational Economics
# Taught by JF Houde

# Alex von Hafften

using DataFrames, StatFiles, LinearAlgebra, Statistics

cd("/Users/alexandervonhafften/Documents/UW Madison/problem_sets/econ_899b/ps3")

# Question 1 - Invert demand for 1985

# define matrices
characteristics = Float64.(Array(select(filter(row -> row.Year == 1985, df_characteristics), :price, :dpm, :hp2wt, :size, :turbo, :trans, :model_class_2, :model_class_3, :model_class_4, :model_class_5, :cyl_2, :cyl_4, :cyl_6, :cyl_8, :drive_2, :drive_3, :Intercept)))
instruments     = Float64.(Array(select(filter(row -> row.Year == 1985, df_instruments), :i_import, :diffiv_ed_0, :diffiv_local_0, :diffiv_local_1, :diffiv_local_2, :diffiv_local_3)))
shares          = Float64.(Array(select(filter(row -> row.Year == 1985, df_characteristics), :share)))
incomes         = Float64.(Array(df_types))
prices          = Float64.(Array(select(filter(row -> row.Year == 1985, df_characteristics), :price)))

# define parameters
λ = 0.6
R = length(incomes)
J = length(prices)

# compute idiosyncratic utility
μ = zeros(R, J)
for i in 1:R
    for j in 1:J
        μ[i, j] = λ * incomes[i] * prices[j]
    end
end


δ = zeros(J)

function demand(δ::Array{Float64}, μ::Array{Float64})
    e = zeros(size(μ))
    for r in 1:size(μ)[1]
        for j in 1:size(μ)[2]
            e[r,j] = exp(δ[j] + μ[r, j])
        end
    end
end

















# IV regression function
function ivreg(mY::Array{Float64}, mVariables::Array{Float64}, mInstruments::Array{Float64}, mWeight::Array{Float64})  
    mInvXZX=inv((mVariables'* mInstruments)*mWeight*(mInstruments'* mVariables))
    if mInvXZX==0
        vIV=fill!(zeros(size(mVariables)[2], 1), NaN)
    else 
        vIV=mInvXZX*(mVariables'mInstruments)*mWeight*(mInstruments'mY)
    end
    return vIV
end

ivreg(shares_1985, characteristics_1985, instruments_1985, Matrix(1.0I, 6, 6))


