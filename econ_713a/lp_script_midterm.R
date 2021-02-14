
# DAA
library(matchingMarkets)

m.prefs <- matrix(c(1,2,3,
                    2,1,3,
                    3,2,1),
                  byrow = FALSE, ncol = 3, nrow = 3)

m.prefs

w.prefs <- matrix(c(3,2,1,
                    1,3,2,
                    3,1,2),
                  byrow = FALSE, ncol = 3, nrow = 3)

w.prefs

iaa(s.prefs = m.prefs,
    c.prefs = w.prefs,
    acceptance = "deferred")

iaa(s.prefs = w.prefs,
    c.prefs = m.prefs,
    acceptance = "deferred")

# Linear programming solver
library(lpSolve)

f.obj <- c(0, 0, 0, 0, -1, 0)

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

f.dir <- c(">=",
           "==",
           ">=",
           "==",
           ">=",
           ">=",
           ">=",
           ">=",
           "==")

f.rhs <- c(3,
           4,
           4,
           7,
           6,
           4,
           5,
           5,
           8)

lp("max", f.obj, f.con, f.dir, f.rhs)$solution