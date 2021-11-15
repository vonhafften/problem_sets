# Computational Economics
# Professor JF Houde
# Problem set 2
# Alex von Hafften 
# November 22, 2021

############################################################
# one-dimensional integration
############################################################

@everywhere using CSV, DataFrames

# quadrature nodes and weights
@everywhere KPU_1d = DataFrame(CSV.File("PS2/KPU_d1_l20.csv"))

@everywhere function integrate_1d(f, upper_bound)

    # define functions to translate the (0, 1) interval into appropriate interval
    points = log.(KPU_1d[:, :x1]) .+ upper_bound
    jacobian = 1 ./KPU_1d[:, :x1]

    # sum over grid points
    return sum(KPU_1d[:, :weight] .* f.(points) .* jacobian)

end

############################################################
# two-dimensional integration
############################################################

# quadrature nodes and weights
@everywhere KPU_2d = DataFrame(CSV.File("PS2/KPU_d2_l20.csv"))

@everywhere function integrate_2d(f, upper_bound_0, upper_bound_1)

    points_0 = log.(KPU_2d[:, :x1]) .+ upper_bound_0
    jacobian_0 = 1 ./ KPU_2d[:, :x1]

    points_1 = log.(KPU_2d[:, :x2]) .+ upper_bound_1
    jacobian_1 = 1 ./KPU_2d[:, :x2]

    return sum(KPU_2d[:, :weight] .* f.(points_0, points_1) .* jacobian_0 .* jacobian_1)

end
