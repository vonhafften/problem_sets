
# Draw shocks
Random.seed!(1234);
ε = rand(Normal(0, 1), T) * σ_0;

# Create data
x_0 = zeros(T);
x_0[1] = ε[1];
for t = 2:T
    x_0[t] = ρ_0 * x_0[t - 1] + ε[t]
end

# compute M_T
M_T = zeros(3);
M_T[1] = sum(x_0) / T
M_T[2] = sum((x_0 .- M_T[1]).^2) / T
M_T[3] = sum((x_0 .- M_T[1]).* (prepend!(x_0[1:T-1], 0) .- M_T[1])) / T
