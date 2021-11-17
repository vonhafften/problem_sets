# Computational Economics
# Professor JF Houde
# Problem set 2
# Alex von Hafften 
# November 22, 2021

############################################################
# Mortgage data
############################################################

using Tables, DataFrames, CSV, StatFiles

cd("/Users/alexandervonhafften/Documents/UW Madison/problem_sets/econ_899b/ps2/")

df = DataFrame(load("PS2/Mortgage_performance_data.dta"))

# Create new columns.
df[!, :c] .= 1.0

# Define x, y, and z.
df_x = select(df, 
              :c, :score_0, :rate_spread, :i_large_loan, :i_medium_loan, 
              :i_refinance, :age_r, :cltv, :dti, :cu,  :first_mort_r, :i_FHA,
              :i_open_year2, :i_open_year3, :i_open_year4, :i_open_year5)
df_t = select(df, :duration)
df_z = select(df, :score_0, :score_1, :score_2)    

# Define matrices.
x = Float64.(Array(df_x))
t = Float64.(Array(df_t))
z = Float64.(Array(df_z))

# Get dimensions.
N = size(x)[1]
K_x = size(x)[2]
K_z = size(z)[2]


