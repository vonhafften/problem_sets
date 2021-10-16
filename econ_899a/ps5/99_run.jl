####################################################################################################################
# Krusell and Smith (1998)
# Alex von Hafften
# October 13, 2021

# ECON 899A Computational Economics
# Problem set 5

# This file contains runs the model.
####################################################################################################################

include("04_solve_model.jl");

@elapsed results = Solve_model()
