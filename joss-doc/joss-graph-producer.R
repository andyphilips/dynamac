# Series simulator to make pretty pictures from dynardl
library(dynamac)

beta.x <- 2 # coef on x
beta.lx <- -3 # coef on l.x
phi <- 0.5 # coef on ldv
const <- 2.5 # constant needed to have a value of mean.y
mean.y <- const/(1-phi) # steady-state of Y

set.seed(100)
x.full <- rnorm(500)
y.full <- rep(NA, length(x.full))

y.full[1] <- const/1-phi

for(i in 2:length(x.full)) {
#	y.full[i] <- const*rnorm(1, 1, 0.05) + phi*y.full[i - 1] + beta.x*x.full[i] + beta.lx*x.full[i - 1] + rnorm(1, 0, 0.05)
  y.full[i] <- const + phi*y.full[i - 1] + beta.x*x.full[i] + beta.lx*x.full[i - 1] + rnorm(1, 0, 2)
  
}

# Trim the burn in
x <- x.full[101:500]
y <- y.full[101:500]

par(mfrow = c(2, 1))
plot(y)
plot(x)
dev.off()

set.seed(1)
model <- dynardl(y ~ x, lags = list("x" = c(1)), levels = c("x"), 
	simulate = TRUE, fullsims = TRUE, 
	shockvar = "x",
	sims = 10000, range = 20,
	ec = FALSE)
	
summary(model)



############ EC
# set.seed(1)
# const <- 2 # constant needed to have a value of mean.y
# x.full <- rep(NA, 500)
# y.diff.full <- rep(NA, length(x.full))
# y.full <- rep(NA, length(x.full))
# y.full[1] <- 1
# x.full[1] <- 15
# y.diff.full[1] <- beta.lx/(-1 * phi)

set.seed(1)
x.error <- rnorm(500, 0, 2)
x.full <- x.error
y.error <- rnorm(500, 0, 2)
y.full <- y.error
y.diff.full <- rep(NA, length(y.full))
phi <- -0.8 # coef on ldv
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

par(mfrow = c(2, 2))
plot(y.diff)
plot(y)
plot(x)
dev.off()

set.seed(1)
model.ec <- dynardl(y ~ x, lags = list("x" = c(1)), diffs = c("x"), 
	simulate = TRUE, fullsims = TRUE, 
	shockvar = "x",
	sims = 10000, range = 20,
	ec = TRUE)
	
summary(model.ec)

dynardl.all.plots(model.ec, bw = TRUE, tol = 0.05)
dev.off()

beta.lx/(-1 * phi)




### Pairing the types together (of the six) for the levels/ec models (png)
# Figures 1 - 6

par(mfrow = c(1, 2), mar = c(4.1, 4.1, 1.1, 1.1))
dynardl.simulation.plot(model, bw = TRUE, response = "levels", main = "Y in Levels")
dynardl.simulation.plot(model.ec, bw = TRUE, response = "levels", main = "Y in Differences")
dev.off()

par(mfrow = c(1, 2), mar = c(4.1, 4.1, 1.1, 1.1))
dynardl.simulation.plot(model, bw = TRUE, response = "diffs", main = "Y in Levels")
dynardl.simulation.plot(model.ec, bw = TRUE, response = "diffs", main = "Y in Differences")
dev.off()

par(mfrow = c(1, 2), mar = c(4.1, 4.1, 1.1, 1.1))
dynardl.simulation.plot(model, bw = TRUE, response = "levels.from.mean", main = "Y in Levels")
dynardl.simulation.plot(model.ec, bw = TRUE, response = "levels.from.mean", main = "Y in Differences")
dev.off()

par(mfrow = c(1, 2), mar = c(4.1, 4.1, 1.1, 1.1))
dynardl.simulation.plot(model, bw = TRUE, response = "cumulative.diffs", main = "Y in Levels")
dynardl.simulation.plot(model.ec, bw = TRUE, response = "cumulative.diffs", main = "Y in Differences")
dev.off()

par(mfrow = c(1, 2), mar = c(4.1, 4.1, 1.1, 1.1))
dynardl.simulation.plot(model, bw = TRUE, response = "cumulative.abs.diffs", main = "Y in Levels", abs.errors = "within.period")
dynardl.simulation.plot(model.ec, bw = TRUE, response = "cumulative.abs.diffs", main = "Y in Differences", abs.errors = "within.period")
dev.off()

par(mfrow = c(1, 2), mar = c(4.1, 4.1, 1.1, 1.1))
dynardl.simulation.plot(model, bw = TRUE, response = "shock.effect.decay", main = "Y in Levels")
dynardl.simulation.plot(model.ec, bw = TRUE, response = "shock.effect.decay", main = "Y in Differences")
dev.off()






















####### REPLICATION
library(readstata13)

data <- read.dta13("/Users/scj0014/Myfiles/Dropbox/Dynpss Project/Dynamic TS Effects and QOI/Illustrative QOI Manuscript text/Plots/PardosSagarzazu_JOP.dta")


quantile(data$econcurrent_agr, c(0.95), na.rm = T) - 
quantile(data$econcurrent_agr, c(0.05), na.rm = T)


quantile(data$media, c(0.95), na.rm = T) - 
quantile(data$media, c(0.05), na.rm = T)


set.seed(09262019)
model.govtattn <- dynardl(gov_saliency ~ unemp + gdp_growth + mip_imputed + econcurrent_agr + opp_mainparty + media + timereal,
	levels = c("unemp", "gdp_growth", "mip_imputed", "econcurrent_agr", "opp_mainparty", "media", "timereal"),
	lags = list("gov_saliency" = 1, "unemp" = 1, "gdp_growth" = 1, "mip_imputed" = 1, 
		"econcurrent_agr" = 1, "opp_mainparty" = 1, "media" = 1),
	data = data, simulate = TRUE, fullsims = TRUE, shockvar = "econcurrent_agr", sims = 10000, shockval = 1.26, burnin = 100 # This is moving from 5% to 95%, illustrative
)


model.oppattn <- dynardl(opp_mainparty ~ unemp + gdp_growth + mip_imputed + econcurrent_agr + gov_saliency + media + timereal,
	levels = c("unemp", "gdp_growth", "mip_imputed", "econcurrent_agr", "gov_saliency", "media", "timereal"),
	lags = list("gov_saliency" = 1, "unemp" = 1, "gdp_growth" = 1, "mip_imputed" = 1, 
		"econcurrent_agr" = 1, "opp_mainparty" = 1, "media" = 1),
	data = data, simulate = TRUE, fullsims = TRUE, shockvar = "media", sims = 10000, shockval = 0.058, burnin = 100
)




# Table 1

govtattn.tab <- data.frame(matrix(rep(NA, length(summary(model.govtattn$model)$coefficients[,1]) * 4), ncol = 4))
for(i in 1:length(summary(model.govtattn$model)$coefficients[,1])) {
	govtattn.tab[((i)*2 - 1), 1] <- rownames(summary(model.govtattn$model)$coefficients)[i]
	govtattn.tab[((i)*2 - 1), 2] <- round(summary(model.govtattn$model)$coefficients[i,1], digits = 3)
	govtattn.tab[((i)*2), 2] <- paste0("(", round(summary(model.govtattn$model)$coefficients[i,2], digits = 3), ")")	
	govtattn.tab[((i)*2 - 1), 3] <- rownames(summary(model.oppattn$model)$coefficients)[i]
	govtattn.tab[((i)*2 - 1), 4] <- round(summary(model.oppattn$model)$coefficients[i,1], digits = 3)
	govtattn.tab[((i)*2), 4] <- paste0("(", round(summary(model.oppattn$model)$coefficients[i,2], digits = 3), ")")	

}

# Figures 7 and 8

par(mar = c(4.1, 4.1, 1.1, 1.1))
dynardl.all.plots(model.govtattn, bw = TRUE, abs.errors = "within.period")
dev.off()


par(mar = c(4.1, 4.1, 1.1, 1.1))
dynardl.all.plots(model.oppattn, bw = TRUE, tol = 0.01, abs.errors = "within.period")
dev.off()
