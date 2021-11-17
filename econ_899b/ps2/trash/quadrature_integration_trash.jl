# Computational Economics
# Professor JF Houde
# Problem set 2
# Alex von Hafften 
# November 22, 2021

############################################################
# one-dimensional integration
############################################################

# quadrature nodes and weights
KPU_1d = DataFrame(CSV.File("PS2/KPU_d1_l20.csv"))

function integrate_1d(f, a = missing, b = missing)

    u = KPU_1d[:, :x1]

    # define functions to translate the (0, 1) interpointsal into appropriate interval
    if !ismissing(a) & !ismissing(b)
        points = (b - a) .* u .+ a
        jacobian = fill(b - a, size(KPU_1d)[1])
    elseif ismissing(a) & ismissing(b)
        points = -log.((1 .-u)./u)
        jacobian = -(u./(1 .-u)) .* ((-1)./u .- (1 .-u)./(u.^2))
    elseif !ismissing(a) & ismissing(b)
        points = -log.(1 .- u) .+ a
        jacobian = 1/(1 .-u)
    elseif ismissing(a) & !ismissing(b)
        points = log.(u) .+ b
        jacobian = 1 ./u
    else
        error("Error in one-dimensional quadrature integration.")
    end

    # sum over grid points
    return sum(KPU_1d[:, :weight] .* f.(points) .* jacobian)

end

function integrate_1d_functions(f; a = missing, b = missing)

    result = 0

    # define functions to translate the (0, 1) interval into appropriate interval
    # plus rho_p, the jacobian of the translation 
    rho(u) = u
    rho_p(u) = 1
        
    if !ismissing(a) & !ismissing(b)
        rho(u) = (b - a) * u + a
        rho_p(u) = b - a
    elseif ismissing(a) & ismissing(b)
        rho(u) = -log((1-u)/u)
        rho_p(u) = -(u/(1-u)) * ((-1)/u - (1-u)/(u^2))
    elseif !ismissing(a) & ismissing(b)
        rho(u) = -log(1 - u) + a
        rho_p(u) = 1/(1-u)
    elseif ismissing(a) & !ismissing(b)
        rho(u) = log(u) + b
        rho_p(u) = 1/u
    else
        error("Error in one-dimensional quadrature integration.")
    end

    # sum over grid points
    for i = 1:size(KPU_1d)[1]
        result += KPU_1d[i, :weight] * f(rho(KPU_1d[i, :x1])) * rho_p(KPU_1d[i, :x1])
    end

    return result

end

############################################################
# two-dimensional integration
############################################################

# quadrature nodes and weights
KPU_2d = DataFrame(CSV.File("PS2/KPU_d2_l20.csv"))

function integrate_2d(f, b_0, b_1)

    result = 0
    
    rho_0(u) = log(u) + b_0
    rho_p_0(u) = 1/u
            
    rho_1(u) = log(u) + b_1
    rho_p_1(u) = 1/u

    for i = 1:size(KPU_2d)[1]
        result += KPU_2d[i, :weight] * f(rho_0(KPU_2d[i, :x1]), rho_1(KPU_2d[i, :x2])) * rho_p_0(KPU_2d[i, :x1]) * rho_p_1(KPU_2d[i, :x2])
    end

    return result

end

function integrate_2d(f; a_0 = missing, b_0 = missing, a_1 = missing, b_1 = missing)
    u_0 = KPU_2d[:, :x1]
    u_1 = KPU_2d[:, :x2]

    if !ismissing(a_0) & !ismissing(b_0)
        points_0 = (b_0 - a_0) .* u_0 .+ a_0
        jacobian_0 = fill(b_0 - a_0, size(KPU_2d)[1])
    elseif ismissing(a_0) & ismissing(b_0)
        points_0 = -log((1 .-u_0)./u_0)
        jacobian_0 = -(u_0./(1 .-u_0)) * ((-1)./ u_0 .- (1 .-u_0) ./(u_0.^2))
    elseif !ismissing(a_0) & ismissing(b_0)
        points_0 = -log(1 .- u_0) + a_0
        jacobian_0 = 1/(1 .-u_0)
    elseif ismissing(a_0) & !ismissing(b_0)
        points_0 = log.(u_0) .+ b_0
        jacobian_0 = 1 ./u_0
    else
        error("Error in two-dimensional quadrature integration.")
    end

    # translation function and jacobian for second dimension
    if !ismissing(a_1) & !ismissing(b_1)
        points_1 = (b_1 - a_1) .* u_1 .+ a_1
        jacobian_1 = fill(b_1 - a_1, size(KPU_2d)[1])
    elseif ismissing(a_1) & ismissing(b_1)
        points_1 = -log((1 .- u_1)./ u_1)
        jacobian_1 = -(u_1 ./(1 .- u_1)) .* ((-1)./u_1 .- (1 .- u_1)./(u_1 .^2))
    elseif !ismissing(a_1) & ismissing(b_1)
        points_1 = -log(1 .- u_1) .+ a_1
        jacobian_1 = 1 ./(1 .-u_1)
    elseif ismissing(a_1) & !ismissing(b_1)
        points_1 = log.(u_1) .+ b_1
        jacobian_1 = 1 ./u_1
    else
        error("Error in two-dimensional quadrature integration.")
    end

    return sum(KPU_2d[:, :weight] .* f.(points_0, points_1) .* jacobian_0 .* jacobian_1)

end
