# Problem Set 3 - BLP
# Computational Economics
# Taught by JF Houde

# Alex von Hafften

# load datasets
df_characteristics = DataFrame(load("PS3/Car_demand_characteristics_spec1.dta"));
df_instruments = DataFrame(load("PS3/Car_demand_iv_spec1.dta"));
df_types = Float64.(DataFrame(load("PS3/Simulated_type_distribution.dta")).Var1);

# Sort datasets
sort!(df_characteristics, [:Year, :Model_id])
sort!(df_instruments, [:Year, :Model_id])

# create X and Z
varlist = ["price","dpm","hp2wt","size","turbo","trans","Year_1986","Year_1987","Year_1988","Year_1989","Year_1990","Year_1991","Year_1992","Year_1993","Year_1994",
		"Year_1995","Year_1996","Year_1997","Year_1998","Year_1999","Year_2000","Year_2001","Year_2002","Year_2003","Year_2004","Year_2005",
		"Year_2006","Year_2007","Year_2008","Year_2009","Year_2010","Year_2011","Year_2012","Year_2013","Year_2014","Year_2015","model_class_2",
		"model_class_3","model_class_4","model_class_5","cyl_2","cyl_4","cyl_6","cyl_8","drive_2","drive_3","Intercept"]

exo_varlist = ["dpm","hp2wt","size","turbo","trans","Year_1986","Year_1987","Year_1988","Year_1989","Year_1990","Year_1991","Year_1992","Year_1993","Year_1994",
		    "Year_1995","Year_1996","Year_1997","Year_1998","Year_1999","Year_2000","Year_2001","Year_2002","Year_2003","Year_2004","Year_2005",
		    "Year_2006","Year_2007","Year_2008","Year_2009","Year_2010","Year_2011","Year_2012","Year_2013","Year_2014","Year_2015","model_class_2",
		    "model_class_3","model_class_4","model_class_5","cyl_2","cyl_4","cyl_6","cyl_8","drive_2","drive_3","Intercept"]

instr_varlist = ["i_import","diffiv_local_0","diffiv_local_1","diffiv_local_2","diffiv_local_3","diffiv_ed_0"]

X = Matrix(df_characteristics[!, varlist])

Z = hcat(Matrix(df_characteristics[!, exo_varlist]), Matrix(df_instruments[!, instr_varlist]))

markets = unique(df_characteristics.Year)