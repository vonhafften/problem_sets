

# Computes μ using T_star operator to verify matrix algebra
function T_operator(R::Results, P::Primitives, G::Grids)
    μ_p = zeros(G.k_n, G.z_n)

    for (i_k, k) in enumerate(G.k_grid), (i_z, z) in enumerate(G.z_grid)
        for (i_k_p, k_p) in enumerate(G.k_grid), (i_z_p, z_p) in enumerate(G.z_grid)
            μ_p[i_k_p, i_z_p] += R.μ[i_k, i_z]  * (R.k_pf[i_k, i_z] == k_p) * R.x_pf[i_k, i_z] * G.z_Π[i_z, i_z_p] 
            μ_p[i_k_p, i_z_p] += R.B * G.ϕ[i_z] * (R.k_pf[1,   i_z] == k_p) * R.x_pf[1,   i_z] * G.z_Π[i_z, i_z_p]
        end
    end

    return μ_p
end


function compute_μ(R::Results)
    P = Primitives()
    G = Grids()

    err, i, maxiter = 100, 1, 100

    while (err > 0.01) & (i < maxiter)
        μ_p = T_operator(R, P, G) 
        err = maximum(abs.(R.μ .- μ_p))
        
        # println("Iteration #", i) # for debugging
        # println("Error is ", err) # for debugging
        
        R.μ = μ_p
        i += 1
    end

    println("μ converged in ", i, " iterations.")

    return R
end