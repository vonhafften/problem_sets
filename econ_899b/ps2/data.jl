# Computational Economics
# Professor JF Houde
# Problem set 2
# Alex von Hafften 
# November 22, 2021

df = DataFrame(load("PS2/Mortgage_performance_data.dta"))

# create new columns
df[!, :c] .= 1.0

# Define y and x
df_y = select(df, :duration)
df_x = select(df, :c, :score_0, :rate_spread, :i_large_loan, :i_medium_loan, 
              :i_refinance, :age_r, :cltv, :dti, :cu,  :first_mort_r, :i_FHA,
              :i_open_year2, :i_open_year3, :i_open_year4, :i_open_year5)
df_z = select(df, :score_0, :score_1, :score_2)    

y = Float64.(Array(df_y))
x = Float64.(Array(df_x))
z = Float64.(Array(df_z))

N = size(x)[1]
K = size(x)[2]