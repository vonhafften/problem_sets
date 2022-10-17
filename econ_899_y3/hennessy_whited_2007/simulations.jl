


using Statistics

# structure for Simulation_Results
mutable struct Simulation_Results

    N::Int64            # number of simulations
    T::Int64            # periods
    burn_in::Int64      # number of period to burn-in simulations
    lz::Matrix{Float64} # log productivity
    w::Matrix{Float64}  # net worth
    k::Matrix{Float64}  # capital
    b::Matrix{Float64}  # debt
    q::Matrix{Float64}  # bond price
    mv::Matrix{Float64} # market value
end

# initialize the simulation results
function Initialize_Simulation_Results()
    N       = 20000
    T       = 14
    burn_in = 200

    lz = zeros(N, burn_in + T + 1)
    w  = zeros(N, burn_in + T + 1)
    k  = zeros(N, burn_in + T + 1)
    b  = zeros(N, burn_in + T + 1)
    q  = zeros(N, burn_in + T + 1)
    mv = zeros(N, burn_in + T + 1)

    return Simulation_Results(N, T, burn_in, lz, w, k, b, q, mv)
end

function simulate_model(P::Primitives, G::Grids, R::Results)
    S = Initialize_Simulation_Results()

    # interpolate model results
    vf_i    = LinearInterpolation((G.grid_w, G.grid_lz), R.vf)
    pf_b_i  = LinearInterpolation((G.grid_w, G.grid_lz), R.pf_b)
    pf_k_i  = LinearInterpolation((G.grid_w, G.grid_lz), R.pf_k)
    q_i     = LinearInterpolation((G.grid_k, G.grid_b, G.grid_lz), R.q)
    w_bar_i = LinearInterpolation(G.grid_lz, R.w_bar)

    # initialize net worth at midpoint
    S.w[:,1] .= (G.min_w + G.max_w)/2

    # iterate over simulated firms
    for i = 1:S.N
        # simulate productivity
        S.lz[i,:] = simulate(G.MC_lz_S, S.burn_in + S.T + 1)
        
        # iterate over time periods
        for t = 1:(S.burn_in + S.T)

            # compute market value
            S.mv[i,t] = vf_i(S.w[i, t], S.lz[i, t])

            # get capital and debt choice from policy functions
            S.b[i, t+1] = pf_b_i(S.w[i, t], S.lz[i, t])
            S.k[i, t+1] = pf_k_i(S.w[i, t], S.lz[i, t])
            S.q[i, t]   = q_i(S.k[i, t+1], S.b[i, t+1], S.lz[i, t])

            # compute taxable income for next period
            y_p = exp(S.lz[i, t+1]) * S.k[i, t+1] ^ P.α - P.δ * S.k[i, t+1] - (1-S.q[i, t]) * S.b[i, t+1]
            
            # compute net worth for next period; max of realized net worth and net worth after defaulting on debt
            S.w[i, t+1] = max(y_p - T_C(y_p, P) + S.k[i, t+1] - S.q[i, t]*S.b[i, t+1], w_bar_i(S.lz[i, t+1]))
        end
    end

    S.w  = S.w[:,(end - S.T):(end-1)]
    S.lz = S.lz[:,(end - S.T):(end-1)]
    S.k  = S.k[:,(end - S.T):(end-1)]
    S.b  = S.b[:,(end - S.T):(end-1)]
    S.q  = S.q[:,(end - S.T ):(end-1)]
    S.mv = S.mv[:,(end - S.T):(end-1)]

    return S
end

function compute_moments(P::Primitives, S::Simulation_Results)
    # investment over book real assets
    inv_ta  = (S.k[:, 2:end] .- (1 .- P.δ) .* S.k[:, 1:(end-1)]) ./ S.k[:, 1:(end-1)]

    # cash flow over book real assets
    y = exp.(S.lz) .* S.k .^ P.α .- P.δ .* S.k .- (1 .- S.q) .* S.b # taxable income
    cf_ta = (exp.(S.lz) .* S.k .^P.α - broadcast(y -> T_C(y, P), y) .+ S.q .* S.b) ./ S.k

    # tobin's q
    tobin_q = (S.mv .+ S.q .* S.b) ./ S.k

    # operating income over book real assets
    oi_ta = exp.(S.lz) .* S.k .^ P.α ./ S.k

    # debt over market value real assets
    d_mv = S.b[:, 2:end] ./ (S.mv[:, 1:(end-1)] .+ S.b[:, 2:end])

    # equity inssurance over book real assets
    ei_ta = max.(-(S.k[:,2:end] .- S.w[:,1:(end-1)] .- S.b[:,2:end]) ./ S.k[:,1:(end-1)], 0.0)

    # cash distribution over book real assets
    div_ta = max.((S.k[:,2:end] .- S.w[:,1:(end-1)] .- S.b[:,2:end]) ./ S.k[:,1:(end-1)], 0.0)

    # create long data for covariances and correlations
    long_data = zeros((S.T-1)*S.N, 10)
    long_data[:,1] = reshape(((1:S.N) *fill(1, S.T-1)')', (S.T-1)*S.N) # firm identicator
    long_data[:,2] = reshape(inv_ta', (S.T-1)*S.N)             # inv_ta 
    long_data[:,3] = reshape(ei_ta', (S.T-1)*S.N)              # equity issuance
    long_data[:,4] = reshape(d_mv', (S.T-1)*S.N)               # leverage
    long_data[:,5] = reshape(oi_ta[:,2:end]', (S.T-1)*S.N)     # current oi_ta
    long_data[:,6] = reshape(oi_ta[:,1:(end-1)]', (S.T-1)*S.N) # lag oi_ta

    # vector for moments
    results = zeros(12)
   
    # moments
    results[1] = mean(ei_ta)        # average equity issuance over assets
    results[2] = var(ei_ta)         # variance of equity issuance over assets
    results[3] = var(inv_ta)        # variance of investment over assets
    results[4] = mean(ei_ta .> 0.0) # frequency of equity issuances
    results[5] = mean(div_ta)       # mean payout ratios # might be wrong
    results[6] = mean(d_mv .< 0.0)  # frequency of negative debt
    results[7] = var(div_ta)        # variance of distributions
    results[8] = mean(d_mv)         # average debt-assets
    results[9] = cov(long_data[:,2], long_data[:,3]) # covariance of inv and equity issuance
    results[10] = cov(long_data[:,2], long_data[:,4]) # covariance of inv and leverage
    results[11] = cor(long_data[:,5], long_data[:,6]) # serial correlation in operating income over total assets
    results[12] = std(long_data[:,5] .- results[11] * long_data[:,6]) # standard deviation of the shock
    
    return results
end