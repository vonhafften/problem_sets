
setwd("/Users/vonhafften/Documents/UW Madison/problem_sets/econ_717b/ps2/")

png("p1_q3_a.png", width = 480, height = 360)
plot(1, type = "n", ylim = c(-2, 2), xlim = c(-2, 2), xlab = "U_V", ylab= "V")

abline(h=0, v=0)

abline(a = 0, b= 1, col = "red")
abline(a = 1, b= 1, col = "blue")

text("Z = 1", x = 1, y = 2, pos = 4, col = "blue")
text("Z = 0", x = -1.5, y = -1.5, pos = 2, col = "red")

abline(v = 0, col = "purple", lty = 2, lwd= 3)
abline(v = -1, col = "forestgreen", lty = 2, lwd = 3)

text("Always Takers", x = 1, y = -1)
text("Compliers", x= -0.5, y = 1)
text("Never Takers", x= -1.5, y = 1)

dev.off()
