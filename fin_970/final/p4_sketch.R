
setwd("/Users/vonhafften/Documents/UW Madison/problem_sets/fin_970/final/")
lambda <- 0.1

N <- rpois(100, lambda)

mu <- -1

Q <- N*mu - lambda*mu

png("p4_1.png", width = 480, height = 360)
plot(Q, type = "l", xlab = "", ylab = "Q")
dev.off()

png("p4_2.png", width = 480, height = 360)
plot(cumsum(Q), type = "l", xlab = "", ylab = "Cumulative Sum of Q")
dev.off()
