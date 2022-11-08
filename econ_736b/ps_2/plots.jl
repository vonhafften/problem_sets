using Plots

grid_p = 0.5:0.02:0.9
grid_q_tilde = 0.9:0.005:0.99
grid_R = 1.1:0.1:3.0
grid_I = 0.0:0.5:10.0
grid_α = 0.2:0.02:0.7

N_p = length(grid_p)
N_q_tilde = length(grid_q_tilde)
N_R = length(grid_R)
N_I = length(grid_I)
N_α = length(grid_α)

j_p = Int64(floor(N_p/2))
j_q_tilde = Int64(floor(N_q_tilde/2))
j_R = Int64(floor(N_R/2))
j_I = Int64(floor(N_I/2))
j_α = Int64(floor(N_α/2))

π_1_results = zeros(N_p, N_q_tilde, N_R, N_α, N_I);
π_2_results = zeros(N_p, N_q_tilde, N_R, N_α, N_I);
π_3_results = zeros(N_p, N_q_tilde, N_R, N_α, N_I);
π_4_results = zeros(N_p, N_q_tilde, N_R, N_α, N_I);
cases = fill(0, N_p, N_q_tilde, N_R, N_α, N_I);

π_1(p, q_tilde, R, α, I) = 0
π_2(p, q_tilde, R, α, I) = p^2 * ((p^2 + q_tilde - 1)^(α/(1-α)) * R^(2/(1-α))*(α^(α/(1-α)) - α^(1/(1-α))) + I/(p^2 + q_tilde - 1))
π_3(p, q_tilde, R, α, I) = R^(1/(1-α)) * (p^2 + 2*p*(1-p) + q_tilde - 1)^(α/(1-α)) * ((p^2 * R + 2*p*(1-p)) * α^(α/(1-α)) * ((p^2 * R + 2*p*(1-p))/(p^2 + 2*p*(1-p))^(α/(1-α))) - (p^2 + 2*p*(1-p))*α^(1/(1-α))*((p^2*R + 2*p*(1-p))/(p^2 + 2*p*(1-p)))^(1/(1-α)) ) + (p^2 + 2*p*(1-p))*I/(p^2+2*p*(1-p) + q_tilde - 1)
π_4(p, q_tilde, R, α, I) = (p^2 * R^2 + 2 * p * (1-p)*R + (1-p)^2)^(1/(1-α)) * q_tilde^(α/(1-α)) * (α^(α/(1-α))- α^(1/(1-α))) + I/q_tilde

for (i_p, p) = enumerate(grid_p), (i_q_tilde, q_tilde) = enumerate(grid_q_tilde),  (i_R, R) = enumerate(grid_R), (i_α, α) = enumerate(grid_α), (i_I, I) = enumerate(grid_I)
    # println("==============")
    # println("p = ", p)
    # println("q_tilde = ", q_tilde)
    # println("R = ", R)
    # println("α = ", α)
    # println("I = ", I)

    π_1_results[i_p, i_q_tilde, i_R, i_α, i_I] = π_1(p, q_tilde, R, α, I)
    π_2_results[i_p, i_q_tilde, i_R, i_α, i_I] = π_2(p, q_tilde, R, α, I)
    π_3_results[i_p, i_q_tilde, i_R, i_α, i_I] = π_3(p, q_tilde, R, α, I)
    π_4_results[i_p, i_q_tilde, i_R, i_α, i_I] = π_4(p, q_tilde, R, α, I)

    max_π = maximum([π_1_results[i_p, i_q_tilde, i_R, i_α, i_I], π_2_results[i_p, i_q_tilde, i_R, i_α, i_I], π_3_results[i_p, i_q_tilde, i_R, i_α, i_I], π_4_results[i_p, i_q_tilde, i_R, i_α, i_I]])
    
    if π_4_results[i_p, i_q_tilde, i_R, i_α, i_I] == max_π
        cases[i_p, i_q_tilde, i_R, i_α, i_I] = 4
    elseif  π_3_results[i_p, i_q_tilde, i_R, i_α, i_I] == max_π
        cases[i_p, i_q_tilde, i_R, i_α, i_I] = 3
    elseif  π_2_results[i_p, i_q_tilde, i_R, i_α, i_I] == max_π
        cases[i_p, i_q_tilde, i_R, i_α, i_I] = 2
    else
        cases[i_p, i_q_tilde, i_R, i_α, i_I] = 1
    end
end 

total = N_p * N_q_tilde* N_R*  N_α*  N_I
sum(cases .== 4)/total
sum(cases .== 3)/total
sum(cases .== 2)/total
sum(cases .== 1)/total

plot(grid_p, π_2_results[:, j_q_tilde, j_R, j_α, j_I], label = "Produce if ε = R²")
plot!(grid_p, π_3_results[:, j_q_tilde, j_R, j_α, j_I], label = "Produce if ε ≥ R")
plot!(grid_p, π_4_results[:, j_q_tilde, j_R, j_α, j_I], label = "Produce if ε ≥ 1")
title!(string("p = ", grid_p[j_p], ", q̃ = ", grid_q_tilde[j_q_tilde], ", R = ", grid_R[j_R], ", α = ", grid_α[j_α], ", I = ", grid_I[j_I]))
xaxis!("p")
yaxis!("π")

plot(grid_q_tilde, π_2_results[j_p, :, j_R, j_α, j_I], label = "Produce if ε = R²")
plot!(grid_q_tilde, π_3_results[j_p, :, j_R, j_α, j_I], label = "Produce if ε ≥ R")
plot!(grid_q_tilde, π_4_results[j_p, :, j_R, j_α, j_I], label = "Produce if ε ≥ 1")
title!(string("p = ", grid_p[j_p], ", q̃ = ", grid_q_tilde[j_q_tilde], ", R = ", grid_R[j_R], ", α = ", grid_α[j_α], ", I = ", grid_I[j_I]))
xaxis!("q̃")
yaxis!("π")

plot(grid_R, π_2_results[j_p, j_q_tilde, :, j_α, j_I], label = "Produce if ε = R²")
plot!(grid_R, π_3_results[j_p, j_q_tilde, :, j_α, j_I], label = "Produce if ε ≥ R")
plot!(grid_R, π_4_results[j_p, j_q_tilde, :, j_α, j_I], label = "Produce if ε ≥ 1")
title!(string("p = ", grid_p[j_p], ", q̃ = ", grid_q_tilde[j_q_tilde], ", R = ", grid_R[j_R], ", α = ", grid_α[j_α], ", I = ", grid_I[j_I]))
xaxis!("q̃")
yaxis!("π")

plot(grid_α, π_2_results[j_p, j_q_tilde, j_R, :, j_I], label = "Produce if ε = R²")
plot!(grid_α, π_3_results[j_p, j_q_tilde, j_R, :, j_I], label = "Produce if ε ≥ R")
plot!(grid_α, π_4_results[j_p, j_q_tilde, j_R, :, j_I], label = "Produce if ε ≥ 1")
title!(string("p = ", grid_p[j_p], ", q̃ = ", grid_q_tilde[j_q_tilde], ", R = ", grid_R[j_R], ", α = ", grid_α[j_α], ", I = ", grid_I[j_I]))
xaxis!("α")
yaxis!("π")

plot(grid_I, π_2_results[j_p, j_q_tilde, j_R, j_α, :], label = "Produce if ε = R²")
plot!(grid_I, π_3_results[j_p, j_q_tilde, j_R, j_α, :], label = "Produce if ε ≥ R")
plot!(grid_I, π_4_results[j_p, j_q_tilde, j_R, j_α, :], label = "Produce if ε ≥ 1")
title!(string("p = ", grid_p[j_p], ", q̃ = ", grid_q_tilde[j_q_tilde], ", R = ", grid_R[j_R], ", α = ", grid_α[j_α], ", I = ", grid_I[j_I]))
xaxis!("I")
yaxis!("π")