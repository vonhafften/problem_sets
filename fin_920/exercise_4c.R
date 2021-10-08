
# choices of a_l and a_h
a_l <- 0:10
b_l <- 0:10

# uses E(a) = E(b) = 5 and Var(a) = 2 and Var(b) = 1 to get p_a and p_b
p_a <- 2/(a_l^2 - 10 * a_l +27)
p_b <- 2/(b_l^2 - 10 * b_l +26)

# uses E(a) = E(b) = 5 and p_a and p_b to get a_h and b_h
a_h <- (5 - a_l * p_a) / (1 - p_a)
b_h <- (5 - b_l * p_b) / (1 - p_b)

# calculates expected utility of a and b
k = 10
u_a <- p_a * log(k + a_l) + (1 - p_a) * log(k + b_h)
u_b <- p_b * log(k + b_l) + (1 - p_b) * log(k + b_h)

plot(y=u_a, x= a_l, type = "l")
lines(y=u_b, x = b_l, col = "red")

u_a > u_b

e_a_3 = p_a * a_l^3 + (1-p_a) * a_h^3
e_b_3 = p_b * b_l^3 + (1-p_b) * b_h^3

e_a_3 > e_b_3

i <- 4

a_l[i]
a_h[i]
p_a[i]
u_a[i]

b_l[i]
b_h[i]
p_b[i]
u_b[i]

plot_range <- c(min(k + a_l[i], k + b_l[i]), max(k + a_h[i], k + b_h[i]))

plot(x = plot_range, y= log(plot_range), type = "l")

abline(v = k + a_l[i], col = "red")
abline(v = k + a_h[i], col = "red")
abline(h = u_a[i], col = "red")

abline(v = k + b_l[i], col = "blue")
abline(v = k + b_h[i], col = "blue")
abline(h = u_b[i], col = "blue")
