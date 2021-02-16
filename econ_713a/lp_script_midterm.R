# Linear programming solver
library(lpSolve)

f.obj <- c(0, 0, -1, 0, 0, 0)

f.con <- matrix(c(1, 0, 0, 1, 0, 0,
                  0, 1, 0, 1, 0, 0,
                  0, 0, 1, 1, 0, 0,
                  1, 0, 0, 0, 1, 0,
                  0, 1, 0, 0, 1, 0,
                  0, 0, 1, 0, 1, 0,
                  1, 0, 0, 0, 0, 1,
                  0, 1, 0, 0, 0, 1,
                  0, 0, 1, 0, 0, 1),
                nrow = 9,
                byrow = TRUE)

f.dir <- c("==",
           ">=",
           ">=",
           ">=",
           "==",
           ">=",
           ">=",
           ">=",
           "==")

f.rhs <- c(10,
           15,
           17,
           0,
           10,
           14,
           0,
           0,
           10)

lp("max", f.obj, f.con, f.dir, f.rhs)$solution
