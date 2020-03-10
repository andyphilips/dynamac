library(MASS)
library(lmtest)
source('R/dynamac.R')


set.seed(1)
x.error <- rnorm(500, 0, 2)
x.full <- x.error
y.error <- rnorm(500, 0, 2)
y.full <- y.error
y.diff.full <- rep(NA, length(y.full))
phi <- 0 #-0.8 # coef on ldv
beta.lx <- 1 # coef on l.x
beta.diff.x <- -2 # coef on x


for(i in 2:length(x.full)) {
  x.full[i] <- x.full[i - 1] + x.error[i]
  y.diff.full[i] <- phi*y.full[i - 1] + beta.diff.x*(x.full[i] - x.full[i - 1]) + beta.lx*x.full[i - 1] + y.error[i] # rnorm(1, 0, 1) + const*rnorm(1, 1, 0.25) +
  y.full[i] <- y.full[i - 1] + y.diff.full[i]
}

# Trim the burn in
x <- x.full[101:500]
y.diff <- y.diff.full[101:500]
y <- y.full[101:500]

# the default: adds lagged dependent variable if specified with ec = TRUE
set.seed(1)
source('R/dynamac.R')

model.ec <- dynardl(y ~ x, lags = list("x" = c(1)), diffs = c("x"), 
                    simulate = TRUE, fullsims = TRUE, 
                    shockvar = "x",
                    sims = 10000, range = 20,
                    ec = TRUE)
summary(model.ec)

par(mfrow = c(1, 4))
plot(y, type = "l")
plot(dshift(y), type = "l")
dynardl.simulation.plot(model.ec, response = "diffs")
dynardl.simulation.plot(model.ec, response = "levels")


# the mod: do not add lagged dependent variable if specified with ec = TRUE, mimicks formula (11) in Philips 2018, American J of Political Science, p. 234
set.seed(1)
source('R/dynamac.R')


model.ec.nlag <- dynardl(y ~ x, lags = list("x" = c(1)), diffs = c("x"), 
                    simulate = T, fullsims = TRUE, 
                    shockvar = "x",
                    sims = 10000, range = 20,
                    ec = TRUE, ec_enforce_lagged_dv = F)
summary(model.ec.nlag)

par(mfrow = c(1, 4))
plot(y, type = "l")
plot(dshift(y), type = "l")

dynardl.simulation.plot(model.ec.nlag, response = "diffs")
dynardl.simulation.plot(model.ec.nlag, response = "levels")

plot(model.ec.nlag$simulation$central+5039.383,type='l')
dynardl.simulation.plot(model.ec, response = "levels")

