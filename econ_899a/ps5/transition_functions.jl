# Krusell and Smith (1998)
# Alex von Hafften
# October 13, 2021

# ECON 899A Computational Economics
# Problem set 5

# Based on matlab available on Dean's teaching website: 
# https://sites.google.com/a/wisc.edu/deancorbae/teaching

# Recession transition function
function get_Π_z()
    # parameters of transition matrix:
    durgd = 8.0
    durbd = 8.0

    Π_z = zeros(2, 2)

    # z transition probabilities
    Π_z[1, 1]   = (durgd-1)/durgd
    Π_z[1, 2]   = 1 - (durbd-1)/durbd
    Π_z[2, 1]   = 1 - (durgd-1)/durgd
    Π_z[2, 2]   = (durbd-1)/durbd
    
    Π_z
end

# The employment transition function 
function get_Π_ε()
    Π_z = get_Π_z()

    durug = 1.5
    unempg= 0.04
    unempb= 0.1
    durub = 2.5

    # ε transition probabilities 
    pgg00 = (durug-1)/durug
    pbb00 = (durub-1)/durub
    pbg00 = 1.25*pbb00  
    pgb00 = 0.75*pgg00
    pgg01 = (unempg - unempg*pgg00)/(1-unempg)
    pbb01 = (unempb - unempb*pbb00)/(1-unempb)
    pbg01 = (unempb - unempg*pbg00)/(1-unempg)
    pgb01 = (unempg - unempb*pgb00)/(1-unempb)
    pgg10 = 1 - (durug-1)/durug
    pbb10 = 1 - (durub-1)/durub
    pbg10 = 1 - 1.25*pbb00  
    pgb10 = 1 - 0.75*pgg00
    pgg11 = 1 - (unempg - unempg*pgg00)/(1-unempg)
    pbb11 = 1 - (unempb - unempb*pbb00)/(1-unempb)
    pbg11 = 1 - (unempb - unempg*pbg00)/(1-unempg) 
    pgb11 = 1 - (unempg - unempb*pgb00)/(1-unempb) 

    # matrix
    Π_ε      = zeros(4, 4)
    Π_ε[1,1] = Π_z[1, 1]*pgg11
    Π_ε[2,1] = Π_z[2, 1]*pbg11
    Π_ε[3,1] = Π_z[1, 1]*pgg01
    Π_ε[4,1] = Π_z[2, 1]*pbg01
    Π_ε[1,2] = Π_z[1, 2]*pgb11
    Π_ε[2,2] = Π_z[2, 2]*pbb11
    Π_ε[3,2] = Π_z[1, 2]*pgb01
    Π_ε[4,2] = Π_z[2, 2]*pbb01
    Π_ε[1,3] = Π_z[1, 1]*pgg10
    Π_ε[2,3] = Π_z[2, 1]*pbg10
    Π_ε[3,3] = Π_z[1, 1]*pgg00
    Π_ε[4,3] = Π_z[2, 1]*pbg00
    Π_ε[1,4] = Π_z[1, 2]*pgb10
    Π_ε[2,4] = Π_z[2, 2]*pbb10
    Π_ε[3,4] = Π_z[1, 2]*pgb00
    Π_ε[4,4] = Π_z[2, 2]*pbb00

    Π_ε
end

# Compute stationary distribution of any Markov process
function compute_Π_star(Π::Array{Float64})    
    x = length(Π[1,:])
    
    Π_star = zeros(x)
    Π_tmp  = Π^10000

    for i=1:x
        Π_star[i] = Π_tmp[i,i]
    end

    Π_star
end
