# Problem Set 3 - BLP
# Computational Economics
# Taught by JF Houde

# Alex von Hafften

using DataFrames, StatFiles, LinearAlgebra, Statistics, Plots

cd("/Users/alexandervonhafften/Documents/UW Madison/problem_sets/econ_899b/ps3")

# read in functions
include("data.jl")
include("model.jl")

####################################################################################
# Question 1 - Invert demand for 1985 and plot sup norm error
####################################################################################

err_list_1985 = invert_demand(1985, 0.6)[2]
err_list_1985_cm = invert_demand(1985, 0.6; tolerence = [1e-12,1e-12])[2]

# plot omitting first iteration with huge error
plot(err_list_1985[2:end], label = "CM + NM");
hline!([1], label="CM + NM Threshold");
plot!(err_list_1985_cm[2:end], label = "CM");
xlabel!("Error");
ylabel!("Iteration #")

savefig("question_1.png")











characteristics = Float64.(Array(select(filter(row -> row.Year == 1985, df_characteristics), :price, :dpm, :hp2wt, :size, :turbo, :trans, :model_class_2, :model_class_3, :model_class_4, :model_class_5, :cyl_2, :cyl_4, :cyl_6, :cyl_8, :drive_2, :drive_3, :Intercept)))
instruments     = Float64.(Array(select(filter(row -> row.Year == 1985, df_instruments), :i_import, :diffiv_ed_0, :diffiv_local_0, :diffiv_local_1, :diffiv_local_2, :diffiv_local_3)))
