
using Parameters, Plots, Statistics

cd("/Users/vonhafften/Documents/UW Madison/problem_sets/econ_810a/project/")

include("model.jl")

P = Primitives()

vf_m_e_i = LinearInterpolation((P.grid_a, P.grid_b), R.vf_m[P.T, 1, :, :])
vf_m_u_i = LinearInterpolation((P.grid_a, P.grid_b), R.vf_m[P.T, 2, :, :])
vf_u_e_i = LinearInterpolation((P.grid_a, P.grid_b), R.vf_u[P.T, 1, :, :])
vf_u_u_i = LinearInterpolation((P.grid_a, P.grid_b), R.vf_u[P.T, 2, :, :])

# figure out price
function price_function(p, a_p, b_p, b_hat)
    if (p < 0)  | (b_p + b_hat > P.max_b)
        return -1
    else
        if (a_p + p + 0.001 >= P.max_a) & (a_p + p < P.max_a)
            epsilon = -0.001
        elseif (a_p + p >= P.max_a)
            return -1
        else 
            epsilon = 0.001
        end

        # get value for intermediary
        W_b_hat = W(b_hat, P.δ, R.β_I)

        # compute value function at a' + p and b'
        V_1 = R.γ * P.Π_y[i_y, 1] * vf_m_e_i[a_p + p, b_p]
        V_1 += R.γ * P.Π_y[i_y, 2] * vf_m_u_i[a_p + p, b_p]
        V_1 += (1 - R.γ) * P.Π_y[i_y, 1] * vf_u_e_i[a_p + p, b_p]
        V_1 += (1 - R.γ) * P.Π_y[i_y, 2] * vf_u_u_i[a_p + p, b_p]

        # compute value function at a' + p + epsilon and b'
        V_1_eps = R.γ * P.Π_y[i_y, 1] * vf_m_e_i[a_p + p + epsilon, b_p]
        V_1_eps += R.γ * P.Π_y[i_y, 2] * vf_m_u_i[a_p + p + epsilon, b_p]
        V_1_eps += (1 - R.γ) * P.Π_y[i_y, 1] * vf_u_e_i[a_p + p + epsilon, b_p]
        V_1_eps += (1 - R.γ) * P.Π_y[i_y, 2] * vf_u_u_i[a_p + p + epsilon, b_p]

        # get numerical derivative at a' + p and b'
        V_1_p = (V_1_eps - V_1)/epsilon

        # compute value function at a' and b' + b_hat
        V_2 = R.γ * P.Π_y[i_y, 1] * vf_m_e_i[a_p, b_p + b_hat]
        V_2 += R.γ * P.Π_y[i_y, 2] * vf_m_u_i[a_p, b_p + b_hat]
        V_2 += (1 - R.γ) * P.Π_y[i_y, 1] * vf_u_e_i[a_p, b_p + b_hat]
        V_2 += (1 - R.γ) * P.Π_y[i_y, 2] * vf_u_u_i[a_p, b_p + b_hat]

        return P.θ * V_1_p * (b_hat - p) - (1-P.θ) * (V_1 - V_2)                        
    end
end



########
R_1 = Initialize(0.99, 0.2)
Solve_terminal_period!(R_1)
prices = 0:0.1:2

i_y = 2

price_conditions_2_5_1 = price_function.(prices, 2.0, 4.0, 1.0)
price_conditions_2_5_2 = price_function.(prices, 2.0, 4.0, 2.0)
price_conditions_2_5_3 = price_function.(prices, 2.0, 4.0, 3.0)

plot(R.vf_u[P.T, 1, [1, 10, P.N_a], :]', label=["low cash" "med cash" "high cash"], legend = :bottomright)
title!("Employed, unmatched, t = 30")
xlabel!("LT bond holdings")
ylabel!("value")

plot(R.vf_u[P.T, 2, [1, 10, P.N_a], :]', label=["low cash" "med cash" "high cash"], legend = :bottomright)
title!("Unemployed, unmatched, t = 30")
xlabel!("LT bond holdings")
ylabel!("value")

plot(prices, price_conditions_2_5_1, label = "a' = 2.0, b' = 5.0, b_hat' = 1.0")
plot!(prices, price_conditions_2_5_2, label = "a' = 2.0, b' = 5.0, b_hat' = 2.0")
plot!(prices, price_conditions_2_5_3, label = "a' = 2.0, b' = 5.0, b_hat' = 3.0")
plot!(prices, zeros(length(prices)), label = false, color = :black)
title!("Pricing Condition (t = 29 and y = 0.5)")
xlabel!("P")
ylabel!("g(P)")
savefig("pricing_condition.png")