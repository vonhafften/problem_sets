#keyword-enabled structure to hold model primitives
@with_kw struct Primitives
    β::Float64 = 0.99 #discount rate
    δ::Float64 = 0.025 #depreciation rate
    α::Float64 = 0.36 #capital share
    k_min::Float64 = 0.01 #capital lower bound
    k_max::Float64 = 75.0 #capital upper bound
    nk::Int64 = 1000 #number of capital grid points
    k_grid::Array{Float64,1} = collect(range(k_min, length = nk, stop = k_max)) #capital grid
    ζ::Array{Float64, 1} = [1.25, 0.2]
    Π::Array{Float64, 2} = [0.977 0.023; 0.074 0.926]
    nz::Int64 = length(ζ)
end

#structure that holds model results
mutable struct Results
    val_func::Array{Float64, 2} #value function
    pol_func::Array{Float64, 2} #policy function
end

#function for initializing model primitives and results
function Initialize()
    prim = Primitives() #initialize primtiives
    val_func = [zeros(prim.nk) zeros(prim.nk)] #initial value function guess
    pol_func = [zeros(prim.nk) zeros(prim.nk)] #initial policy function guess
    res = Results(val_func, pol_func) #initialize results struct
    prim, res #return deliverables
end

#Bellman Operator
function Bellman(prim::Primitives, res::Results)
    @unpack val_func = res #unpack value function
    @unpack k_grid, β, δ, α, nk, ζ, Π, nz = prim #unpack model primitives
    v_next = [zeros(nk) zeros(nk)] #next guess of value function to fill

    for k_index = 1:nk
        k = k_grid[k_index] #value of k

        for z_i = 1:nz

            candidate_max = -Inf #bad candidate max
            budget = ζ[z_i] * k^α + (1 - δ) * k #budget

            for kp_index in 1:nk #loop over possible selections of k'
                c = budget - k_grid[kp_index] #consumption given k' selection
                if c > 0 #check for positivity
                    val = log(c)  #compute value
                    for zp_i = 1:nz
                        val = val + β * Π[z_i, zp_i] * val_func[kp_index, zp_i]
                    end
                    if val > candidate_max #check for new max value
                        candidate_max = val #update max value
                        res.pol_func[k_index, z_i] = k_grid[kp_index] #update policy function
                    end
                end
            end
            v_next[k_index, z_i] = candidate_max #update value function
        end
    end
    v_next #return next guess of value function
end

#Value function iteration
function V_iterate(prim::Primitives, res::Results; tol::Float64 = 1e-4, err::Float64 = 100.0)
    n = 0 #counter

    while err>tol #begin iteration
        v_next = Bellman(prim, res) #spit out new vectors
        err = abs.(maximum(v_next.-res.val_func))/abs(v_next[prim.nk, 1]) #reset error level
        res.val_func = v_next #update value function
        n+=1
    end
    println("Value function converged in ", n, " iterations.")
end

#solve the model
function Solve_model(prim::Primitives, res::Results)
    V_iterate(prim, res) #in this case, all we have to do is the value function iteration!
end
##############################################################################
