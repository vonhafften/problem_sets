# Alex von Hafften 
# Besanko and Doraszelski
# Supplement to Computational Economics
# December 23, 2021

# This file contains the plotting functions

############################################################################
# Plots
############################################################################

using Plots

function plot_surface(matrix::Array{Float64})
    P = Primitives()
    function f(x, y)
        matrix[x, y]
    end
    plot(1:P.n_q, 1:P.n_q, f, st=:surface, legend = false)
    xlabel!("q_bar_1")
    ylabel!("q_bar_2")
end

function plot_static_results(R::Results)

    p1 = plot_surface(R.q_star_1);
    title!("q_star_1");

    p2 = plot_surface(R.q_star_2);
    title!("q_star_2");

    p3 = plot_surface(R.p_star_1);
    title!("p_star_1");

    p4 = plot_surface(R.p_star_2);
    title!("p_star_2");

    p5 = plot_surface(R.π_1);
    title!("π_1");

    p6 = plot_surface(R.π_2);
    title!("π_2");

    plot(p1, p2, p3, p4, p5, p6, layout =  (3, 2));
    plot!(size=(700,900), titlefontsize = 12, guidefontsize=8)
end

function plot_dynamic_results(R::Results)

    p1 = plot_surface(R.x_pf_1);
    title!("x_pf_1");

    p2 = plot_surface(R.x_pf_2);
    title!("x_pf_2");

    p3 = plot_surface(R.vf_1);
    title!("vf_1");

    p4 = plot_surface(R.vf_2);
    title!("vf_2");

    plot(p1, p2, p3, p4, layout =  (2, 2));
    plot!(size=(700,700), titlefontsize = 12, guidefontsize=8)
end