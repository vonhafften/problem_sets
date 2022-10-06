# Code to solve Huggett model in cts time 
# a la Achdou et al (2017)

# Alex von Hafften
# Sept. 22

using Parameters, LinearAlgebra, SparseArrays, Arpack

# structure for model primitives
@with_kw struct Primitives

    # Model parameters
    λ::Float64   = 0.05  # earning shock arrival rate
    γ::Float64   = 2.5   # coefficient of relative risk aversion
    ρ::Float64   = 0.018 # discrount rate

    # grids 
    # earnings
    y_L::Float64            = 0.01       # low earning
    y_H::Float64            = 0.03       # high earning
    y_grid::Vector{Float64} = [y_L, y_H] # save as vector
    N_y::Int64              = 2          # number of earning states

    # other
    Δ_im::Float64 = 1000.0  # implicit method convergence step
end

# structure to hold results
mutable struct Results
    vf::Matrix{Float64}    # value function
    pf::Matrix{Float64}    # consumption policy function
    a_dot::Matrix{Float64} # asset drift implied by consumption policy function

    dvf::Matrix{Float64}    # value function derivative
    dvf_F::Matrix{Float64}  # value function forward derivative
    dvf_B::Matrix{Float64}  # value function backward derivative
    dvf_0::Matrix{Float64}  # value function stay-put derivative
    I_F::Matrix{Float64}    # indicator for value function forward derivative
    I_B::Matrix{Float64}    # indicator for value function backward derivative
    I_0::Matrix{Float64}    # indicator for value function stay-put derivative

    π::Matrix{Float64}        # stationary distribution
    r::Float64                # interest rate
    net_supply_bonds::Float64 # net supply of bonds

    # assets
    min_a::Float64              # minimum
    max_a::Float64              # maximum
    N_a::Int64                  # number of grid points
    a_grid_srl::Vector{Float64} # grid as step range length
    a_grid::Vector{Float64}     # grid
    Δ_a::Float64                # step size
end

# initialize Results structure
function Initialize(r::Float64, bc::Float64)
    P = Primitives()

    # asset grid
    min_a      = bc  
    max_a      = 5.0      
    N_a        = 1001   
    a_grid_srl = range(min_a, max_a; length=N_a)
    a_grid     = collect(a_grid_srl)
    Δ_a        = a_grid[2]-a_grid[1]

    vf    = zeros(N_a, P.N_y)
    pf    = zeros(N_a, P.N_y)
    a_dot = zeros(N_a, P.N_y)
    dvf   = zeros(N_a, P.N_y)
    dvf_F = zeros(N_a, P.N_y)
    dvf_B = zeros(N_a, P.N_y)
    dvf_0 = zeros(N_a, P.N_y)
    I_F   = zeros(N_a, P.N_y)
    I_B   = zeros(N_a, P.N_y)   
    I_0   = zeros(N_a, P.N_y)   
    π     = zeros(N_a, P.N_y)    
    net_supply_bonds = 0.0


    Results(vf, pf, a_dot, dvf, dvf_F, dvf_B, dvf_0, I_F, I_B, I_0, π, r, net_supply_bonds, min_a, max_a, N_a, a_grid_srl, a_grid, Δ_a)
end

# utility
function u(c::Float64, γ::Float64)
    if c > 0.0
        return (c^(1 - γ))/(1 - γ)
    else
        -1/eps()
    end
end

# marginal utility
function mu(c::Float64, γ::Float64)
    if c > 0.0
        return c^(- γ)
    else 
        1/eps()
    end
end

# marginal utility inverse
function mu_inv(mu::Float64, γ::Float64)
    mu^(-1/γ)
end

# solve model
function Solve!(R::Results)
    P = Primitives()

    # fill in initial values
    for (i_a, a) = enumerate(R.a_grid), (i_y, y) = enumerate(P.y_grid)
        # initial guess for value function is staying put forever
        R.vf[i_a, i_y] = u(a*R.r + y, P.γ)/P.ρ

        # stay put derivative doesn't change
        R.dvf_0[i_a, i_y] = mu(a*R.r + y, P.γ)
    end

    # solve value function
    error = 100
    n = 1
    tol = 1e-4
    while error > tol
        upwind(R)
        error = update_vf(R)

        #println("Iteration #", n)
        #println("Current error: ", error)
        n += 1
    end

    #println("Completed in iteration: ", n, ".")
    #println("Final error: ", error)

    compute_π(R)
end

# upwind procedure
function upwind(R::Results)
    P = Primitives()

    # reset indicator functions
    R.I_F = zeros(R.N_a, P.N_y)
    R.I_B = zeros(R.N_a, P.N_y)
    R.I_0 = zeros(R.N_a, P.N_y)

    # loop over state space
    for (i_a, a) = enumerate(R.a_grid), (i_y, y) = enumerate(P.y_grid)

        # compute forward derivative
        if i_a == R.N_a
            R.dvf_F[i_a, i_y] = mu(y + R.r*a, P.γ)
        else
            R.dvf_F[i_a, i_y] = (R.vf[i_a+1, i_y] - R.vf[i_a, i_y])/R.Δ_a
        end

        # compute backward derivative
        if i_a == 1
            R.dvf_B[1, i_y] = mu(y + R.r*a, P.γ)
        else
            R.dvf_B[i_a, i_y] = (R.vf[i_a, i_y] - R.vf[i_a-1, i_y])/R.Δ_a
        end

        # compute a_dot based on forward and backward derivatives
        c_F     = mu_inv(R.dvf_F[i_a, i_y], P.γ)
        c_B     = mu_inv(R.dvf_B[i_a, i_y], P.γ)
        a_dot_F = y + R.r*a - c_F
        a_dot_B = y + R.r*a - c_B

        # figure out indicators
        if a_dot_F > 0
            R.I_F[i_a, i_y] = 1.0
        elseif a_dot_B < 0
            R.I_B[i_a, i_y] = 1.0
        elseif (a_dot_F < 0) & (a_dot_B > 0)
            R.I_0[i_a, i_y] = 1.0
        end
    end

    # hard-coding backward change at top end of asset grid
    R.I_B[R.N_a,:] .= 1.0
    R.I_F[R.N_a,:] .= 0.0
    R.I_0[R.N_a,:] .= 0.0

    # update derivative
    R.dvf = R.I_0 .* R.dvf_0 + R.I_F .* R.dvf_F + R.I_B .* R.dvf_B
    
    # compute policy functions 
    # consumption
    R.pf  = mu_inv.(R.dvf, P.γ)
    # implies change in assets
    for (i_a, a) = enumerate(R.a_grid), (i_y, y) = enumerate(P.y_grid)
        R.a_dot[i_a, i_y] = R.r * a + y - R.pf[i_a, i_y]
    end
end

# function for updating value function using implicit method
function update_vf(R::Results)
    P = Primitives()

    # initalize next value function iteration
    vf_next = zeros(R.N_a, P.N_y)
    
    # compute A - transition matrix
    A = compute_A(R)

    # compute B matrix
    B = (P.ρ + (1/P.Δ_im)) * Matrix{Float64}(I, P.N_y*R.N_a, P.N_y*R.N_a) .- A

    # compute b matrix
    b = u.(R.pf, P.γ) + (1/P.Δ_im) .* R.vf
    b_stack = reshape(b, (R.N_a*P.N_y, 1))
    
    # vf_next 
    vf_next_stacked = B\b_stack
    vf_next = reshape(vf_next_stacked, (R.N_a, P.N_y))

    # sup norm
    error = maximum(abs.(vf_next - R.vf))

    # update value function
    R.vf = vf_next
    
    return error
end

function compute_π(R::Results)
    P = Primitives()

    # compute A
    A = compute_A(R)
    AT = A'

    #doctor some entires to avoid singularity
    temp_col = zeros(R.N_a*P.N_y)
    temp_col[1]=0.1
    temp_row = zeros(1,R.N_a*P.N_y)
    temp_row[1] = 1
    AT[1,:] = temp_row

    #compute the distribution
    π_stack = (A'\temp_col)
    π_mass = ones(1,P.N_y*R.N_a) * π_stack * R.Δ_a
    π_stack = π_stack/π_mass
    R.π = reshape(π_stack, (R.N_a, P.N_y))

    #compute the net supply of bonds
    net_supply = 0.0
    for i_y = P.N_y
        net_supply -= sum(R.π[:, i_y] .* R.a_grid) * R.Δ_a
    end
    R.net_supply_bonds = net_supply
end

function compute_A(R::Results)
    P = Primitives()

    A = spzeros(P.N_y*R.N_a, P.N_y*R.N_a)
    
    for (i_a, a) = enumerate(R.a_grid), (i_y, y) = enumerate(P.y_grid)
        a_dot = R.a_dot[i_a, i_y]

        x = (-R.I_B[i_a, i_y] .* a_dot) ./ R.Δ_a
        y = (R.I_B[i_a, i_y] .* a_dot) ./ R.Δ_a .- (R.I_F[i_a, i_y] .* a_dot) ./ R.Δ_a .- P.λ
        z = R.I_F[i_a, i_y] .* a_dot ./ R.Δ_a
        
        # find no shock index
        j = (i_y-1)*R.N_a + i_a
        A[j, j] = y
        if i_a > 1
            A[j, j-1] = x
        end
        if i_a < R.N_a
            A[j, j+1] = z
        end

        # add lambda is on the diagonal of the off-diagonal blocks for shocks
        for (i_y_not, y_not) = enumerate(P.y_grid)
            if i_y != i_y_not
                k = R.N_a*(i_y_not-1) + i_a
                A[j, k] = P.λ
            end
        end
    end

    return A
end

function Solve!(r_min::Float64, r_max::Float64, bc::Float64)
    R = Initialize(r_min, bc)
    try Solve!(R)
        if R.net_supply_bonds < 0.0
            error("Net supply of bonds is negative at r_min")
        end
    catch
        error("Can't solve problem at r_min")
    end

    R.r = r_max
    try Solve!(R)
        if R.net_supply_bonds > 0.0
            error("Net supply of bonds is positive at r_max")
        end
    catch
        error("Can't solve problem at r_max")
    end


    tol = 0.00001
    while true
        R.r = (r_min + r_max)/2
        println("Interest rate: ", R.r)
        Solve!(R)
        println("Net supply of bonds: ", R.net_supply_bonds)

        if R.net_supply_bonds > tol
            r_min = (r_min + r_max)/2
        elseif  R.net_supply_bonds < -tol
            r_max = (r_min + r_max)/2
        else
            break
        end

    end

    R
end
