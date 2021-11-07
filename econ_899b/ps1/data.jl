# Computational Economics
# Professor JF Houde
# Problem set 1
# Alex von Hafften 
# November 7, 2021

df = DataFrame(load("Mortgage_performance_data.dta"))

# create new columns
df[!, :c] .= 1.0

# Define y and x
df_y = select(df, :i_close_first_year)
df_x = select(df, :c, :i_large_loan, :i_medium_loan, :rate_spread,
              :i_refinance, :age_r, :cltv, :dti, :cu,  :first_mort_r,
              :score_0, :score_1, :i_FHA, :i_open_year2, :i_open_year3,
              :i_open_year4, :i_open_year5)
              
y = Float64.(Array(df_y))
x = Float64.(Array(df_x))

N = size(x)[1]
K = size(x)[2]