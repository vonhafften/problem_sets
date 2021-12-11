# Problem Set 3
# Computational Economics
# Taught by JF Houde

# Alex von Hafften

using DataFrames, StatFiles, LinearAlgebra, Statistics

cd("/Users/alexandervonhafften/Documents/UW Madison/problem_sets/econ_899b/ps3")

include("blp_func_ps.jl")

# Main dataset and panel structure
df_PanelCharact = DataFrame(load("PS3/Car_demand_characteristics_spec1.dta"))
n = size(df_PanelCharact)[1]
vYear = unique(df_PanelCharact.Year)
T = length(vYear)

# Outcome variables
vShare = df_PanelCharact.share
vDelta_iia = df_PanelCharact.delta_iia
vDelta0 = vDelta_iia
vPrice = df_PanelCharact.price

# characteristics and instrumental variables

# Load characteristics
varlist = ["price","dpm","hp2wt","size","turbo","trans","Year_1986","Year_1987","Year_1988","Year_1989","Year_1990","Year_1991","Year_1992","Year_1993","Year_1994",
		"Year_1995","Year_1996","Year_1997","Year_1998","Year_1999","Year_2000","Year_2001","Year_2002","Year_2003","Year_2004","Year_2005",
		"Year_2006","Year_2007","Year_2008","Year_2009","Year_2010","Year_2011","Year_2012","Year_2013","Year_2014","Year_2015","model_class_2",
		"model_class_3","model_class_4","model_class_5","cyl_2","cyl_4","cyl_6","cyl_8","drive_2","drive_3","Intercept"]
exo_varlist = ["dpm","hp2wt","size","turbo","trans","Year_1986","Year_1987","Year_1988","Year_1989","Year_1990","Year_1991","Year_1992","Year_1993","Year_1994",
		    "Year_1995","Year_1996","Year_1997","Year_1998","Year_1999","Year_2000","Year_2001","Year_2002","Year_2003","Year_2004","Year_2005",
		    "Year_2006","Year_2007","Year_2008","Year_2009","Year_2010","Year_2011","Year_2012","Year_2013","Year_2014","Year_2015","model_class_2",
		    "model_class_3","model_class_4","model_class_5","cyl_2","cyl_4","cyl_6","cyl_8","drive_2","drive_3","Intercept"]
mX = Float64.(Array(select(df_PanelCharact, varlist)))
mX_summary = hcat(varlist, mean(mX, dims =1)') # Mean product characteristics

# Load price and differentations IVs
df_DemandIV = DataFrame(load("PS3/Car_demand_iv_spec1.dta"))
aIVlist = ["i_import","diffiv_local_0","diffiv_local_1","diffiv_local_2","diffiv_local_3","diffiv_ed_0"]
mExclIV = Float64.(Array(select(df_DemandIV, aIVlist)))
mExclIV_summary = hcat(aIVlist, mean(mExclIV, dims =1)') # Mean cost IV (import) and differentiation measures.
mIV = hcat(Float64.(Array(select(df_PanelCharact, exo_varlist))), mExclIV)

# Non-linear attributes
mZ = Float64.(Array(select(df_PanelCharact, :price)))

# Pre-compute the row IDs for each market
aProductID = [Vector{Int64}() for i in 1:T];
for i = 1:T
    aProductID[i] = findall(==(vYear[i]),df_PanelCharact.Year)
end

# Random coefficients
mEta = Float64.(DataFrame(load("PS3/Simulated_type_distribution.dta")).Var1)
Sim = length(mEta)

println("Distribution of eta")
println("Mean: ", mean(mEta))
println("SD: ", std(mEta))
println("Minimum: ", minimum(mEta))
println("1st Quartile: ", quantile(mEta, 0.25))
println("Median: ", quantile(mEta, 0.5))
println("3rd Quartile: ", quantile(mEta, 0.75))
println("Maximum: ", maximum(mEta))

# Pre-compute interaction between price and random-coefficient
# Two dimenional arrays of JxSim matrices: T x Nb of variables
aZ = Array{Array{Float64,2},2}(undef, T, size(mZ)[2])
for i = 1:T
    for j = 1:size(mZ)[2]
        aZ[i,j] = mZ[aProductID[i]] .* mEta[:,j]'
    end
end

#################################################################
# GMM Estimator
#################################################################

vParam = zeros(size(mZ)[2])
vParam0 = vParam
step = 0

# 2SLS weighting matrix
A = inv(mIV' * mIV)

println("Plot the iteration process")

# Inversion algorithm
vParam[1] = 0.6
# Contraction mapping
inverse(vDelta0, vParam, 0, 10^(-12), vDelta0, aProductID, vShare)


#### Progress.....




