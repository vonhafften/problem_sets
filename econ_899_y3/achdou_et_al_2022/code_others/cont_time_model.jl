#struct to hold model primitives
@with_kw struct Primitives
    lambda::Float64 = 0.05 #shock arrival rate
    gamma::Float64 = 2.5 #CRRA coefficient
    rho::Float64 = 0.018 #discount rate
    a_min::Float64 = -0.5 #minimum assets
    a_max::Float64 = 5.0 #maximum assets
    Na::Int64 = 1000 #number of asset grid points
    a_grid::Array{Float64,1} = collect(range(a_min, length=Na, stop=a_max)) #asset grid
    Da::Float64 = a_grid[2]-a_grid[1] #step size of a_grid; used for calculating derivatives
    y_grid::Array{Float64, 1} = [0.01, 0.03] #possible states
    Ny = length(y_grid) #number of states
    delta::Float64 = 1000.0 #object for convergence-checking
end

#struct to hold model results
mutable struct Results
    val_func::Array{Float64,2}  #value function
    dV::Array{Float64,2}        #value function derivatives
    adot::Array{Float64,2}      #asset drift terms
    pol_func::Array{Float64,2}  #policy function
    Ib::Array{Float64,2}        #indicators for back derivatives; used for weights
    If::Array{Float64,2}        #indicators for forward derivatives; used for weights
    I0::Array{Float64,2}        #indicators for zero derivatives; used for weights
    dVf::Array{Float64,2}       #forward derivative values
    dVb::Array{Float64,2}       #backward derivative values
    dV0 ::Array{Float64,2}      #zero derivative values
    stat_dist::Array{Float64,2}
    r::Float64                  #interest rate
end

#function to call to solve the model
function Solve_model(; delta::Float64 = 1000.0, imp::Int64 = 0, r::Float64 = 0.01)
    prim = Primitives(delta = delta) #initialize primitives
    @unpack Na, Ny, a_grid, gamma, rho, y_grid = prim

    pol_func = zeros(Na, Ny) #policy function guess
    dV = zeros(Na, Ny) #value function derivatives
    adot = zeros(Na, Ny) #drift terms
    Ib = zeros(Na, Ny) #indicators
    If = zeros(Na, Ny) #indicators
    I0 = zeros(Na, Ny) #indicators
    dVf = zeros(Na, Ny) #drift terms
    dVb = zeros(Na, Ny) #drift terms
    dV0 = zeros(Na, Ny) #drift terms
    stat_dist = zeros(Na, Ny) #stationary distribution
    r = r #interest rate
    val_func = zeros(Na, Ny)
    for i = 1:Na, j = 1:Ny #first guess of val func: staying put forever
        val_func[i,j] = ((a_grid[i]*r + y_grid[j])^(1-gamma))/(rho*(1-gamma))
    end

    #preallocate results vector
    res = Results(val_func, dV, adot, pol_func, Ib, If, I0, dVf, dVb, dV0, stat_dist, r)

    V_iterate(prim, res, imp) #value function iteration
    Compute_statdist(prim,res) #compute stationary distribution
    res #return results
end

#value function iteration function
function V_iterate(prim::Primitives, res::Results, imp::Int64; tol::Float64 = 1e-4)
    error = 100
    n=0

    while error > tol
        Upwind(prim, res) #update derivatives guess
        error = Compute_error(prim, res, imp) #update value function and error term
        n+=1
    end
    println("Completed in iteration: ", n, ". Current error: ", error)
end

#upwinding function
function Upwind(prim::Primitives, res::Results)
    @unpack y_grid, a_grid, Na, Ny, gamma, Da =  prim
    @unpack r, val_func = res
    res.Ib = zeros(Na, Ny) #reset these puppies
    res.If = zeros(Na, Ny)
    res.I0 = zeros(Na, Ny)

    for i = 1:Na, j = 1:Ny #loop over state space
        res.dV0[i,j] = (y_grid[j] + r*a_grid[i])^(-gamma) #compute zero derivative

        try #compute forward derivative
            res.dVf[i,j] = (val_func[i+1, j] - val_func[i, j])/Da
        catch err
            if isa(err, BoundsError) #catch endpoint case
                res.dVf[i,j] = (y_grid[j] + r*a_grid[Na])^(-gamma)
            end
        end

        try #compute backward derivative
            res.dVb[i,j] = (val_func[i, j] - val_func[i-1, j])/Da
        catch err
            if isa(err, BoundsError) #catch endpoint case
                res.dVb[i,j] = (y_grid[j] + r*a_grid[1])^(-gamma)
            end
        end

        cf = res.dVf[i,j]^(-1/gamma) #forward consumption
        cb = res.dVb[i,j]^(-1/gamma) #backward consumption
        af = y_grid[j] + r*a_grid[i] - cf #forward drift
        ab = y_grid[j] + r*a_grid[i] - cb #backward drift

        #apply upwinding to update our results vectors
        if af>0
            res.If[i,j] = 1.0 #indicator
        elseif ab<0
            res.Ib[i,j]=1.0
        elseif af<0 && ab>0 #tiebreak rule!
            res.I0[i,j] = 1.0
        end
    end
    res.Ib[Na,:] = [1.0 1.0] #hard-coding backward drift at the maximal asset value; if you don't do this there will be hell to pay.
    res.If[Na,:] = [0.0 0.0]
    res.I0[Na,:] = [0.0 0.0]

    res.dV = res.If.*res.dVf + res.Ib.*res.dVb + res.I0.*res.dV0 #update results derivatives vector
    res.pol_func = res.dV.^(-1/gamma) #update policy function
    for i = 1:Na, j = 1:Ny
        res.adot[i,j] = r * a_grid[i] - res.pol_func[i,j] + y_grid[j] #update policy function
    end
end

#function for updating the value function using either implicit or explicit method
function Compute_error(prim::Primitives, res::Results, imp::Int64)
    @unpack y_grid, a_grid, Na, Ny, gamma, Da, rho, delta, lambda = prim
    @unpack r, val_func, dV, adot, pol_func, Ib, If, I0 = res
    v_next = zeros(Na, Ny)

    if imp==0 #time for the explicit method!
        for i = 1:Na, j = 1:Ny #loop over state space
            jnot = 1
            if j == 1
                jnot=2
            end

            rhs = ((pol_func[i,j])^(1-gamma))/(1-gamma) + adot[i,j]*dV[i,j] +
            lambda*(val_func[i, jnot] - val_func[i,j]) #right-hand side term
            v_next[i, j] = delta*(rhs-rho*val_func[i,j])+val_func[i,j] #update guess
        end
    elseif imp==1 #time for the implicit method!
        A_mat = Compute_A_mat(prim, res) #get a matrix
        B_mat = Matrix{Float64}(I, Ny*Na, Ny*Na).*(rho + (1/delta)) .- A_mat #B matrix
        b_mat = ((pol_func).^(1-gamma))./(1-gamma) + (1/delta).*val_func #little b matrix
        b_stack = [b_mat[:,1]; b_mat[:,2]] #now stacked; julia does like calling the inv() function on B_mat
        v_next_stacked = B_mat\b_stack #stacked version of the next value function
        v_next = [v_next_stacked[1:Na] v_next_stacked[Na+1:2*Na]] #reshuffle to something more familiar
    end
    error = maximum(abs.(v_next - res.val_func)) #update error
    res.val_func = v_next #update value function
    error #return error
end

#function for computing A matrix; used for implicit method and finding stationary distribution
function Compute_A_mat(prim, res)
    @unpack y_grid, a_grid, Na, Ny, gamma, Da, rho, delta, lambda = prim
    @unpack r, val_func, dV, adot, pol_func, Ib, If, I0 = res
    A_mat = spzeros(Ny*Na, Ny*Na) #initialize big sparse matrix
    lambda_vec = ones(Na).*lambda #lambda vector

    for j = 1:Ny #loop over states
        x_vec = (-Ib[:,j].*adot[:,j])./Da
        y_vec = (Ib[:,j].*adot[:,j])./Da .- (If[:,j].*adot[:,j])./Da .- lambda
        z_vec = (If[:,j].*adot[:,j])./Da

        #get mismatch index for lambda diagonals; not particularly robust
        jnot = 1
        if j == 1
            jnot=2
        end

        #fill in diagonals
        A_mat[(Na*(j-1)+1):(Na*(j-1)+Na),(Na*(j-1)+1):(Na*(j-1)+Na)] = spdiagm(0=>y_vec[:]) + spdiagm(-1=>x_vec[2:Na]) + spdiagm(1 => z_vec[1:Na-1])
        A_mat[(Na*(j-1)+1):(Na*(j-1)+Na),(Na*(jnot-1)+1):(Na*(jnot-1)+Na)] = spdiagm(0=>lambda_vec[:]) #lambdas
    end
    A_mat
end

#function for computing stationary distribution
function Compute_statdist(prim::Primitives, res::Results)
    @unpack Na, Ny, Da = prim
    A_mat = Compute_A_mat(prim, res) #get the A matrix
    AT = A_mat'

    #doctor some entires to avoid singularity
    temp_col = zeros(Na*Ny)
    temp_col[1]=0.1
    temp_row = zeros(1,Na*Ny)
    temp_row[1] = 1
    AT[1,:] = temp_row

    #compute the distribution
    gstack = (AT\temp_col)
    gmass = ones(1,Ny*Na)*gstack*Da;
    gstack = gstack/gmass;
    res.stat_dist = [gstack[1:Na] gstack[Na+1:2*Na]] #reshuffle to get the distribution
end
###########
