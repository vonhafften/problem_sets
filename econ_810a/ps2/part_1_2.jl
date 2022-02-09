# Alex von Hafften
# ECON 810: Advanced Macro
# Carter Braxton

# Problem set 2

using DataFrames, FixedEffectModels, RegressionTables, Plots

cd("/Users/alexandervonhafften/Documents/UW Madison/problem_sets/econ_810a/ps2/")

####################################################
# parameters
####################################################

# dims
T = 11
N_l = 500
N_s = 500
# idiosyncratic component
σ_ε = 1000.0
μ_ε = 0.0
# time component
γ_1 = 30000.0
γ_t = 1000.0
# effect of job loss
β_6 = -9000.0


####################################################
# simulate earning data
# columns are time
# rows are individuals
####################################################

# idiosyncratic earning component
ε = reshape(μ_ε .+ σ_ε .* randn(T*(N_l+N_s)), (N_l + N_s, T))

# time component of earnings
γ = (γ_1 .* ones(N_l + N_s, T)) .+ (γ_t .* ones(N_l + N_s) * transpose(0:(T-1)))

# affect of job loss
d = vcat(hcat(zeros(N_l, 5), ones(N_l, T-5)), zeros(N_s, T))
β = β_6 .* d

# total simulated earning panel
y = ε + γ + β

####################################################
# run distributed lag regression
####################################################

# get time and individual fes
individual = (1:(N_l + N_s)) * transpose(fill(1, T))
time = fill(1, N_l + N_s) * transpose(1:T)

# create distributed lags
d_m4 = vcat(hcat(zeros(N_l, 1), ones(N_l), zeros(N_l, T-1-1)), zeros(N_s, T))
d_m3 = vcat(hcat(zeros(N_l, 2), ones(N_l), zeros(N_l, T-1-2)), zeros(N_s, T))
d_m2 = vcat(hcat(zeros(N_l, 3), ones(N_l), zeros(N_l, T-1-3)), zeros(N_s, T))
d_m1 = vcat(hcat(zeros(N_l, 4), ones(N_l), zeros(N_l, T-1-4)), zeros(N_s, T))
d_0 = vcat(hcat(zeros(N_l, 5), ones(N_l), zeros(N_l, T-1-5)), zeros(N_s, T))
d_p1 = vcat(hcat(zeros(N_l, 6), ones(N_l), zeros(N_l, T-1-6)), zeros(N_s, T))
d_p2 = vcat(hcat(zeros(N_l, 7), ones(N_l), zeros(N_l, T-1-7)), zeros(N_s, T))
d_p3 = vcat(hcat(zeros(N_l, 8), ones(N_l), zeros(N_l, T-1-8)), zeros(N_s, T))
d_p4 = vcat(hcat(zeros(N_l, 9), ones(N_l), zeros(N_l, T-1-9)), zeros(N_s, T))
d_p5 = vcat(hcat(zeros(N_l, 10), ones(N_l)), zeros(N_s, T))

# pull together into dataframe
data = DataFrame(earnings = reshape(y, T*(N_l + N_s)), 
                 individual = reshape(individual, T*(N_l + N_s)), 
                 time = reshape(time, T*(N_l + N_s)), 
                 dm4 = reshape(d_m4, T*(N_l + N_s)), 
                 dm3 = reshape(d_m3, T*(N_l + N_s)), 
                 dm2 = reshape(d_m2, T*(N_l + N_s)), 
                 dm1 = reshape(d_m1, T*(N_l + N_s)), 
                 d0 = reshape(d_0, T*(N_l + N_s)), 
                 dp1 = reshape(d_p1, T*(N_l + N_s)), 
                 dp2 = reshape(d_p2, T*(N_l + N_s)), 
                 dp3 = reshape(d_p3, T*(N_l + N_s)), 
                 dp4 = reshape(d_p4, T*(N_l + N_s)), 
                 dp5 = reshape(d_p5, T*(N_l + N_s)))

# run regression
model = reg(data, @formula(earnings ~ dm4 + dm3 + dm2 + dm1 + d0 + dp1 + dp2 + dp3 + dp4 + dp5 + fe(individual) + fe(time)), save=true)

regtable(model; renderSettings = latexOutput("part_1_2_table.tex"), fixedeffects= ["fe_time_1"])

# beta plot
plot((-4):5, model.coef, legend=false);
xlabel!("k");
title!("β_k");
savefig("part_1_2_beta.png")

# gamma plot
plot(unique(model.fe.time), unique(model.fe.fe_time), legend=false);
title!("γ_t");
xlabel!("t");
savefig("part_1_2_gamma.png")
