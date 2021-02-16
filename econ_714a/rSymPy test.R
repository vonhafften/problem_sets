library("rSymPy")

sympy("x, y, z, alpha, beta = symbols('x y z alpha beta', real = True)")
sympy("i,j,k = symbols('i j k', integer = True, positive = True)")






#######

library(rSymPy)

alpha_ <- Var("alpha")
delta_ <- Var("delta")
beta_ <- Var("beta")

l_t_ <- Var("l_t")
l_t1_ <- Var("l_t1")
i_t_ <- Var("i_t")
k_t_ <- Var("k_t")
k_t1_ <- Var("k_t1")
c_t_ <- Var("c_t")
c_t1_ <- Var("c_t1")

i_bar_ <- Var("i_bar")
k_bar_ <- Var("k_bar")
l_bar_ <- Var("l_bar")
c_bar_ <- Var("c_bar")
g_bar_ <- Var("g_bar")

a_t_ <- Var("a_t")
a_t1_ <- Var("a_t1")
g_t_<- Var("g_t")
g_t1_<- Var("g_t1")
tau_hat_l_t_ <- Var("tau_hat_l_t")
tau_hat_l_t1_ <- Var("tau_hat_l_t1")
tau_hat_i_t_ <- Var("tau_hat_i_t")
tau_hat_i_t1_ <- Var("tau_hat_i_t1")

rho_a_ <- Var("rho_a")
rho_g_ <- Var("rho_g")
rho_l_ <- Var("rho_l")
rho_i_ <- Var("rho_i")

###########################
x_ <- var("x")

sympy("ss1 = c_bar + delta*k_bar - x * (k_bar**alpha) *(l_bar **(1-alpha)) ")

sympy("ss2 = c_bar + (1-alpha) * k_bar ** alpha * l_bar ** (-1-alpha)")

sympy("ss3 = 1 - beta *(alpha * k_bar **(alpha - 1)* l_bar ** (1-alpha) + (1- delta))")

sympy("solve([ss1,ss2, ss3], [c_bar, k_bar, l_bar])")

###########################


# consumption-labor equation
sympy("cl_eq = l_t - alpha/(1+alpha) * k_t - (-1)/(1+alpha) * c_t - 1/(1+alpha) * a_t - (-1)/(1+alpha)*tau_hat_l_t")
sympy("cl_eq1 = l_t1 - alpha/(1+alpha) * k_t1 - (-1)/(1+alpha) * c_t1 - 1/(1+alpha) * a_t1 - (-1)/(1+alpha) * tau_hat_l_t1")

# output equation
sympy("o_eq = i_t - (a_t + alpha * k_t + (1-alpha)*l_t - c_bar*c_t - g_bar*g_t) / i_bar")

# lom equation
sympy("lom_eq = k_t1 - (1-delta) *k_t - delta * i_t")

# euler equation
sympy("e_eq = c_t1 - c_t + tau_hat_i_t - beta * alpha * k_bar **(alpha - 1) * l_bar **(1 - alpha) * (a_t1 + (1-alpha)*(l_t1 - k_t1)) + (1-delta)*tau_hat_i_t1")

# AR-1 processes
sympy("ar_a = a_t1 - rho_a * a_t")
sympy("ar_g = g_t1 - rho_g * g_t")
sympy("ar_l = tau_hat_l_t1 - rho_l * tau_hat_l_t")
sympy("ar_a = tau_hat_i_t1 - rho_i * tau_hat_i_t1")

sympy("solve((cl_eq,cl_eq1, o_eq, lom_eq, e_eq, ar_a, ar_g, ar_l, ar_a), (c_t1, k_t1, c_t, k_t, l_t))")