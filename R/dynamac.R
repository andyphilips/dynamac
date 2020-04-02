# version 0.1.10.9001
# 4/2/2020
# Authors: Soren Jordan, Andrew Q. Philips

# Corrections since previous version:
#	Change to pssbounds to allow, but carefully, different cases of models (through restriction = TRUE/FALSE)
#   Override duplicate input in pssbounds (user + dynardl)
#   Move pssbounds estimation (fstat) to pssbounds function
#	Custom ylim to simulation plots
#	Change default in dynardl to c() for levels and diffs (for clarity)
#	A lot of new unit tests for functions, especially expected warnings

# TO DO: 
#	Do user-specified significances need to be in the plots?
#	Unit testing of plots through vdiffr
#   Future release: permanent shifts as opposed to shocks? Forecasting?
#   Add more autocorrelation tests
#	Automatic lag selection?
#	ts data?

## Deprecated functions file
#' @title Deprecated functions in package \pkg{dynamac}
#' @description The functions listed below are deprecated and will be defunct
#' in the near future. When possible, functions with similar functionality are mentioned.
#' Help pages for deprecated functions are available at \code{help-("deprecated")}.
#' @name dynamac-deprecated
#' @keywords internal
NULL


# Datasets exported: 
#' Data on public concern about economic inequality
#'
#' A dataset from: Wright, Graham. 2017. "The political
#' implications of American concerns about economic inequality."
#' Political Behavior 40(2): 321-346. 
#'
#' @format A data frame with 49 rows and 9 variables:
#' \describe{
#'   \item{year}{Year}
#'   \item{mood}{Public mood liberalism}
#'   \item{urate}{Unemployment rate}
#'   \item{concern}{Concern about economic inequality}
#'   \item{demcontrol}{Democratic control of congress}
#'   \item{incshare10}{Proportion of income of top 10 percent}
#'   \item{csentiment}{Consumer sentiment}
#'   \item{incshare01}{Proportion of income of top 1 percent}
#' }
#' @source \url{http://dx.doi.org/10.7910/DVN/UYUU9G}
#' @docType data
#' @keywords datasets
#' @usage data(ineq)
#' @name ineq
NULL


# Datasets exported: 
#' Data on French Energy Consumption and GDP
#'
#' Data on GDP are from World Bank World Development Indicators. Data
#' on energy consumption are from the PB Statistical Review of World
#' Energy (June 2018).
#'
#' @format A data frame with 53 rows and 4 variables:
#' \describe{
#'   \item{country}{Country}
#'   \item{year}{Year}
#'   \item{lnGDP_cons2010USD}{ln(GDP), constant 2010 US dollars}
#'   \item{lnenergy}{ln(energy consumption), millions tons oil equivalent}
#' }
#' @docType data
#' @keywords datasets
#' @usage data(france.data)
#' @name france.data
NULL

# Datasets exported: 
#' Data on US Supreme Court Approval
#'
#' A dataset from: Durr, Robert H., Andrew D. Martin, and Christina Wolbrecht. 2000. "Ideological divergence and public support for the Supreme Court." American Journal of Policial Science 44(4): 768-776.
#'
#' @format A data frame with 42 rows and 9 variables:
#' \describe{
#'   \item{dcalc}{Supreme Court support}
#'   \item{l_dcalc}{Lagged Supreme Court spport}
#'   \item{iddiv}{Ideological divergence}
#'   \item{mooddev}{Mean deviation of Mood}
#'   \item{dirdev}{Mean deviation of percent liberal decisions}
#'   \item{sg}{Rulings against Solicitor General's amicus briefs}
#'   \item{laws}{Laws declared unconstitutional}
#'   \item{presapp}{Approval of president}
#'   \item{congapp}{Approval of Congress}
#' }
#' @source \url{http://dx.doi.org/10.2307/2669280}
#' @docType data
#' @keywords datasets
#' @usage data(supreme.sup)
#' @name supreme.sup
NULL

## Functions:
## Dependencies: 	MASS (for multivariate normal draws)
#					lmtest (for autocorrelation tests)
#
## Functions included:
# (1) lshift()
#	x = [no default]						vector to lag
#	l = [no default]						number of lags
#
# (2) dshift()
#	x = [no default]						vector to difference
#
# (3) ldshift()
#	x = [no default]						vector to lag difference
# 	l = [no default]						number of lagged differences
#
# (4) dynardl()
#	formula = [no default]				model formula
#	data = [list()]						list of variables
#	lags = [list()]						list of lags, like lags = list("x" = c(1,2) ...)
# 	diffs = [c()]						vector of first diffs, like diffs = c("x", "z" ...)
#	lagdiffs = [list()]					list of variables to lagdiff, like lagdiffs = list("x" = c(1) ...)
# 	levels = [c()]						vector of variables in levels, like c("x", "z" ...)
#	ec = [FALSE] 						should model be error correcting? (i.e. differenced y or not)
#	trend = [FALSE]						include a linear trend 
#	constant = [TRUE]					include a constant
#	noLDV = [FALSE]						do not force a LDV in EC models. highly unrecommended
#	modelout = [FALSE]					print model summary	
#	simulate = [FALSE]					simulate model, or just estimate?
#	shockvar = [list()]					variable to be shocked, like "x"
#	shockval = [sd(data[[shockvar]])]	amount to shock variable
#	time = [10]							time to shock variable
#	qoi = ["mean"]						summarize the response with the mean or the median?
#	forceset = [NULL]					list of variables to be forced to a certain value 
#	range = [20] 						range of simulation
#	burnin = [20]						number of time periods to throw away before beginning reported simulation
#	sims = [1000]						number of simulations to run
#	sig = [95] 							significance for simulations
#	expectedval = [FALSE]				expected value of simulation (averaged errors) or predicted value 
#	fullsims = [FALSE]					store the full, raw simulations in the model object (necessary for some plotting)
#
# (5) area.simulation.plot()				<!!!! DEPRECATED !!!!>
#	x = [no default] 					dynardl object containing a simulation and parameters
#	response = [levels]					should the plot be in levels of Y (levels) or changes from mean of Y (mean.changes)?	
#	bw = [FALSE]						should the plot be in black and white?				
#
# (6) spike.simulation.plot()			<!!!! DEPRECATED !!!!>
#	x = [no default] 					dynardl object containing a simulation and parameters
#	response = [levels]					should the plot be in levels of Y (levels) or changes from mean of Y (mean.changes)?	
#	bw = [FALSE]						should the plot be in black and white?
#
# (7) pssbounds()
#	data = [list()]						a dynardl model object
#	obs = [NULL]						the number of observations in the model time series
#	fstat = [NULL]						the f statistic on the joint first lags of variables
#	tstat = [NULL]						t statistic on LDV in EC model
#	case = [NULL]						case of regression: see documentation
#	k = [NULL]							number of k regressors in first lags
#	digits = 3							digits to report of statistics
#	object.out = [FALSE]					do you want to print all of this into an object?
#
# (8) dynardl.auto.correlated()
#	x = [no default]					a dynardl model object
#	bg.type = ["Chisq"]		 			character string of Chisq or F for the type of bg test used
#	digits = [3]						where to round AIC/BIC
#	order = [NULL]						the order of autocorrelation to test in bgtest
#	object.out = [FALSE]					do you want to print all of this into an object?
#
# (9) summary.dynardl()
#	x = [no default]					a dynardl model object
#
# (10) dynardl.simulation.plot()
#	x = [no default]					a dynardl model object
#	type = ["area"]						should the plot be an area plot or a spike plot?
#	response = ["levels"]				track the response of y in "levels", 
#											"levels.from.mean", period over period "diffs", "cumulative.diffs", 
#											"cumulative.abs.diffs" (like a LRE), or the "shock.effect.decay"
#	last.period = [NULL]				let dynardl decide when to stop calculating absolute cumulative diffs, or truncate
# 	tol = [sd(x$model$y[,1], na.rm = T) * 0.001] tolerance: when dynardl decides to stop calculating abs cumulative diffs
#	bw = [FALSE]						black and white or color?
# 	y.lab = [""]						user-defined y-label or reasonable default
# 	x.lab = [""]						user-defined x-label or reasonable default
#	...								arguments to plot


##########################################
# ------------(1) lshift ----------------#
##########################################
#' Take lag transformation of a series
#' @param x a series to be lagged
#' @param l the number of lags
#' @return the lagged series
#' @details
#' \code{lshift} assumes that the series are ordered, that there is no missing data, and that the time intervals are even
#' @importFrom utils head
#' @author Soren Jordan and Andrew Q. Philips
#' @keywords utilities
#' @examples
#' x.var <- runif(50)
#' l.1.x.var <- lshift(x.var, 1)
#' l.2.x.var <- lshift(x.var, 2)
#' head(x.var)
#' head(l.1.x.var)
#' head(l.2.x.var)
#' @export

lshift <- function(x, l){
	stopifnot(is.numeric(l))
	stopifnot(is.numeric(x))
	out <- c(rep(NA,l), head(x,-l))
	out
}

##########################################
# ------------(2) dshift ----------------#
##########################################
#' Take first difference of a series
#' @param x a series to be differenced
#' @return the differenced series
#' @details
#' \code{dshift} assumes that the series are ordered, that there is no missing data, and that the time intervals are even
#' @author Soren Jordan and Andrew Q. Philips
#' @keywords utilities
#' @examples
#' x.var <- seq(0, 50, 5)
#' d.x.var <- dshift(x.var)
#' head(x.var)
#' head(d.x.var)
#' @export

dshift <- function(x){
	stopifnot(is.numeric(x))
	out <- c(NA, diff(x))
	out
}

##########################################
# ------------(3) ldshift ---------------#
##########################################
#' Take the lagged first difference of a series
#' @param x a series to be differenced
#' @param l the number of lags
#' @return the lagged differenced series
#' @details
#' \code{ldshift} assumes that the series are ordered, that there is no missing data, and that the time intervals are even
#' @importFrom utils head
#' @author Soren Jordan and Andrew Q. Philips
#' @keywords utilities
#' @examples
#' x.var <- runif(50)
#' ld.1.x.var <- ldshift(x.var, 1)
#' ld.2.x.var <- ldshift(x.var, 2)
#' head(x.var)
#' head(ld.1.x.var)
#' head(ld.2.x.var)
#' @export

ldshift <- function(x, l){
	stopifnot(is.numeric(l))
	stopifnot(is.numeric(x))
	out <- c(rep(NA,l), head(c(NA,diff(x)),-l))
	out
}

##########################################
# ------------(4) dynardl ---------------#
##########################################
#' Estimate and simulate ARDL model
#' @description
#' Estimate autoregressive distributed lag models and simulate interesting values (if desired)
#' @param formula a symbolic description of the model to be estimated. ARDL
#' models are estimated using linear regression
#' @param data an optional data frame or list containing the the variables in the model
#' @param lags a list of variables and their corresponding lags to be estimated
#' @param diffs a vector of variables to be differenced. Only first differences are supported
#' @param lagdiffs a list of variables to be included in lagged differences
#' @param levels a vector of variables to be included in levels
#' @param ec estimate model in error-correction form, (i.e., \code{y} appears in first-differences). By default, \code{ec} is set to \code{FALSE}, meaning \code{y} will appear in levels.
#' @param trend include a linear time trend. The default is \code{FALSE}
#' @param constant include a constant. The default is \code{TRUE}
#' @param noLDV do not add a lagged dependent variable (LDV) to ARDL models when omitted in formula (special thanks to Hannes Datta). This is not recommended
#' @param modelout print the regression estimates in the console
#' @param simulate simulate the reponse. Otherwise, just the regression model will be estimated. If \code{simulate = FALSE}, options \code{shockvar}, \code{shockval},  \code{time}, \code{qoi}, \code{forceset}, \code{range}, \code{burnin}, \code{sims}, \code{sig}, \code{expectedval}, and \code{fullsims} are ignored. The default is \code{FALSE} so that users can build models without needing to simulate the results each time. When \code{simulate = TRUE}, users are highly encouraged to set a seed before simulation, as with any stochastic exercise
#' @param shockvar the variable to be shocked in the counterfactual simulation. There is no default
#' @param shockval the amount by which the \code{shockvar} should be shocked. The default is one standard deviation of the shocked variable
#' @param time the time period in the simulation for the variable to be shocked
#' @param qoi summarize the response of the dependent variable with the mean or the median. Although the default is \code{mean}, if there is underlying skew in the distribution, it might be better summarized by \code{median}
#' @param forceset by default, in the simulations, variables in levels will be set to their means; variables in differences will be set to 0. Alternatively, users can set any variable in the model to a different value using a list in \code{forceset}. These values can be any user-defined value, including means, medians, percentiles, or other values of interest 
#' @param range the range of the simulation to be conducted
#' @param burnin the number of time periods to disregard before recording the values. These do not include the \code{range}; in other words, they take place before the \code{range} specified above. Users can increase the number of \code{burnin} periods, but probably should not decrease them. The default is 20
#' @param sims the number of simulations to use in creating the quantities of interest (the response of the dependent variable). The default is 1000
#' @param sig the significance level (1 - \code{p}) that the user wants for the simulations. The default level is 95\% significance (\code{sig = 95})
#' @param expectedval if this is \code{TRUE}, the simulation will record the expected values of across the \code{sims} by averaging errors. The default is \code{FALSE}, since expected values do not account for stochastic error present in the model itself
#' @param fullsims whether all of the raw simulations should be stored in the model object. These are required for some of the more advanced plotting functions, especially those that use the simulations to derive confidence intervals about the size of the period-over-period differences. The default is \code{FALSE}
#' @return \code{dynardl} should always return an estimated model. It may or may not be simulated, according to the user. But the relevant regression output, model residuals (which can be tested for autocorrelation), and simulated response (if created) are stored in a list if the model is assigned to an object
#' @details
#' Estimate an auto-regressive distributed lag model. Moreover, enable a graphical interpretation of the results (through \code{\link{dynardl.simulation.plot}}) by simulating the response of the dependent variable to shocks in one of the regressors, and enable the Pesaran, Shin, and Smith (2001) test for cointegration for error-correction models (through \code{\link{pssbounds}})
#' @importFrom graphics lines plot points polygon segments
#' @importFrom stats AIC BIC as.formula coef lm logLik quantile rchisq rnorm sd sigma vcov median
#' @importFrom utils head flush.console head setTxtProgressBar txtProgressBar
#' @importFrom MASS mvrnorm
#' @author Soren Jordan and Andrew Q. Philips
#' @keywords simulation estimation ardl
#' @examples
#' # Using the inequality data from dynamac
#' ardl.model <- dynardl(concern ~ incshare10 + urate, data = ineq, 
#'        lags = list("concern" = 1, "incshare10" = 1),
#'        diffs = c("incshare10", "urate"), 
#'        ec = TRUE, simulate = FALSE)
#' summary(ardl.model)
#' 
#' # Adding a lagged difference of the dependent variable
#' ardl.model.2 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
#'        lags = list("concern" = 1, "incshare10" = 1),
#'        diffs = c("incshare10", "urate"), 
#'        lagdiffs = list("concern" = 1),
#'        ec = TRUE, simulate = FALSE)
#' summary(ardl.model.2)
#'
#' # Does not work: levels and diffs must appear as a vector
#' \donttest{
#' ardl.model.3 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
#'        lags = list("concern" = 1, "incshare10" = 1),
#'        levels = list("urate" = 1),
#'        diffs = list("incshare10" = 1, "urate" = 1), 
#'        lagdiffs = list("concern" = 1),
#'        ec = TRUE, simulate = FALSE)
#' }
#'
#' ardl.model.3 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
#'        lags = list("concern" = 1, "incshare10" = 1),
#'        levels = c("urate"),
#'        diffs = c("incshare10", "urate"), 
#'        lagdiffs = list("concern" = 1),
#'        ec = TRUE, simulate = FALSE)
#' @export

dynardl <- function(formula,
								data = list(),
								lags = list(),
								diffs = c(),
								lagdiffs = list(),
								levels = c(),
								ec = FALSE,
								trend = FALSE,
								constant = TRUE,
								modelout = FALSE,
								noLDV = FALSE,
								simulate = FALSE,
								shockvar = list(),
								shockval = sd(data[[shockvar]], na.rm = T),
								time = 10,
								qoi = "mean",	
								forceset = NULL,		
								range = 20,
								burnin = 20,
								sims = 1000,
								sig = 95,
								expectedval = FALSE,
								fullsims = FALSE) {
	#################
	# Data controls #
	#################
	# If data not passed, create data frame
	if(identical(data, list())) {
		temp.data.mat <- NULL
		for(i in 1:length(all.vars(formula))) {
			assign(paste(all.vars(formula)[i]), get(all.vars(formula)[i]))
			temp.data.mat <- cbind(temp.data.mat, get(all.vars(formula)[i]))
		}
		colnames(temp.data.mat) <- all.vars(formula)
		data <- data.frame(temp.data.mat)
	}	
	#######################
	# Difference controls #
	#######################
	if(length(diffs)) {
		if(isTRUE(is.list(diffs))) {
			stop("Declare differences as a vector of variables to be first-differenced (i.e. diffs = c('X1', ...))")
		}
		else {
	#	if(!is.list(diffs)) {	# Stop for non-list diffs
	#		stop("'diffs' must be a list")
	#	}
	#	else {
			for(i in 1:length(diffs)) {
	#			if(!(is.numeric(diffs[[i]]))) { # Stop for diffs without diff level
	#				stop("Declare differences with numeric list: diffs = list(\"X\" = c(1))")
	#			}
	#			else {			
	#				for(j in 1:diffs[[i]]) {
	#					if(diffs[[i]][j] > 1) {	# Stop for second or more diffs
	#						stop("Second differences ('diffs' > 1) not supported")
	#					}
	#				}
					# Warn if variables aren't in formula
					if(!(diffs[i] %in% all.vars(formula))) {
							warning(paste(paste("Variable"), paste(diffs[i]), paste("in the diffs list not found in model formula. Not included in differences."), sep = " "))
					}
			}
		}
	#		}
	#	}
	}
	################
	# Lag controls #
	################
	if(length(lags)) {
		if (!is.list(lags)) {
			stop("'lags' must be a list")
		}
		else {
			if(as.character(formula[[2]]) %in% names(lags)) { # Stop for first lag of LDV not in lags
				first.lag.ldv <- FALSE
				for(i in 1:length(names(lags))) {
					if(names(lags)[i] == as.character(formula[[2]])) {
						for(o in 1:length(lags[[i]])) {
							# For the lags of LDV, if it includes the first, we're good
							ifelse(lags[[i]][o] == 1, first.lag.ldv <- TRUE, first.lag.ldv <- first.lag.ldv)
						}
					}
				}
				if(first.lag.ldv == FALSE) {
					if(noLDV == TRUE) {
						warning(paste("ARDL model with no LDV specified. Only finite dynamics are possible through lags of X (which might affect simulations). Are you sure you want this?"))
					} else {
						stop("'lags' must include first lag of dependent variable (LDV) (or specify noLDV = TRUE)")
					}
				}
			}
			# Warn if variables aren't in formula
			for(i in 1:length(names(lags))) {
				if(!(names(lags)[i] %in% all.vars(formula))) {
					warning(paste(paste("Variable"), paste(names(lags)[i]), paste("in the lags list not found in model formula. Not included in lagged levels."), sep = " "))
				}
			}
		}
	}
	##############################
	# Lagged difference controls #
	##############################			
	if(length(lagdiffs)){
		if(!is.list(lagdiffs)) {
			stop("'lagdiffs' must be a list")
		}
		else {
			# Warn if variables aren't in formula
			for(i in 1:length(names(lagdiffs))) {
				if(!(names(lagdiffs)[i] %in% all.vars(formula))) {
					warning(paste(paste("Variable"), paste(names(lagdiffs)[i]), paste("in the lagdiffs list not found in model formula. Not included in lagged differences."), sep = " "))
				}
			}
		}
	}
	###################
	# Levels controls #
	###################
	if(length(levels)){
		if(isTRUE(is.list(levels))) {
			stop("Declare levels as a vector of variables to be included in levels (i.e. levels = c('X1', ...))")
		}
		if(as.character(formula[[2]]) %in% levels) {
			stop("Dependent variable cannot appear in contemporaneous levels in ARDL model.")
		}
		else {
			# Warn if variables aren't in formula
			for(i in 1:length(levels)) {
				if(!(levels[i] %in% all.vars(formula))) {
					warning(paste(paste("Variable"), paste(levels[i]), paste("in the levels list not found in model formula. Not included in levels."), sep = " "))
				}
			}
		}
	}
	#####################
	# Forceset controls #
	#####################	
	if(length(forceset)){
		if(as.character(formula[[2]]) %in% names(forceset)) {
			stop("LDV cannot be forceset in dynamic simulation.")
		}
		else {
			# Warn if variables aren't in formula
			for(i in 1:length(names(forceset))) {
				if(!(names(forceset)[i] %in% all.vars(formula))) {
					warning(paste(paste("Variable"), paste(names(forceset)[i]), paste("in the forceset list not found in model formula. Variable not forced."), sep = " "))
				}
			}
		}
	}
	#######################
	# Simulation controls #
	#######################
	shock.message <- NULL
		#####################
		# Shockvar controls #
		#####################
		# The first set of warnings are actually for all shockvars. Even if a simulation isn't estimated, we need to check the
		#  syntax is correct *providing the user actually has a shockvar*
		if(simulate == FALSE) {
			if(!(identical(deparse(substitute(shockvar)), "list()"))) {	# If the user specified a shockvar of any sort
				if(!(grepl("\"", deparse(substitute(shockvar))))) {	# And it's not in quotations
					stop("'shockvar' must be specified in quotation marks.")		# Because it will break when it searches for shockvars later
				}			
				else {	# Will make lots of other conditions unhappy later
					if(length(shockvar) > 1) {
						stop("You specified more than one shockvar. Only specify one shockvar in a dynamic simulation.")
					}
					if(!(shockvar %in% all.vars(formula))) { # If it's not in the formula, that's a problem
						warning("Your shockvar is not in the model formula, but since you are not simulating, I am ignoring it. Specify modeled variables for dynamic simulations.")
					}
					if(!(shockvar %in% diffs)) { #if its not in diffs
						if(!(shockvar %in% names(lagdiffs))) {
							if(!(shockvar %in% levels)) {
								if(!(shockvar %in% names(lags))) {
									warning("Your shockvar is not found in lags, lagdiffs, levels, or differences, but since you are not simulating, I am ignoring it. Specify modeled variables for dynamic simulations.")
								}
							}
						}
					}		
				}
			}
		}
		else { # if simulate == TRUE
			if(identical(deparse(substitute(shockvar)), "list()")) {	# If the user did NOT specify a shockvar of any sort
				stop("'shockvar' must be specified for a dynamic simulation")
			}
			else {				
				if(!(grepl("\"", deparse(substitute(shockvar))))) {	# And it's not in quotations
					stop("'shockvar' must be specified in quotation marks.")		# Because it will break when it searches for shockvars later
				}			
				else { 
					if(length(shockvar) > 1) {
						stop("Only specify one shockvar in a dynamic simulation")
					}
					else {
						if(!(shockvar %in% all.vars(formula))) { # If it's not in the formula, that's a problem
							stop(paste(paste("Variable"), paste(shockvar), paste("(the shockvar) not found in model formula. Shock variable required for dynamic simulation."), sep = " "))
						}
						else { # If it's not estimated, that's another problem
							if(!(shockvar %in% diffs)) { #if its not in diffs
								if(!(shockvar %in% names(lagdiffs))) {
									if(!(shockvar %in% levels)) {
										if(!(shockvar %in% names(lags))) {
											stop(paste(paste("Variable"), paste(shockvar), paste("(the shockvar) not found lags, lagdiffs, levels, or differences. Shock variable required for dynamic simulation."), sep = " "))
										}
									}
								}
							}
						}
					}	
				}
			}
		if(shockval == sd(data[[shockvar]], na.rm = T)) { # If it's the default shock value, remind the user
			shock.message <- paste(paste(shockvar), paste("shocked by one standard deviation of"), paste(shockvar), paste("by default."), sep = " ")
		}
		#####################
		# Forceset controls #
		#####################	
		if(length(forceset)) {
			if(as.character(formula[[2]]) %in% names(forceset)) {
				stop("LDV cannot be forceset in dynamic simulation.")
			}
			else {
				# Warn if variables aren't in formula
				for(i in 1:length(names(forceset))) {
					if(!(names(forceset)[i] %in% all.vars(formula))) {
						warning(paste(paste("Variable"), paste(names(forceset)[i]), paste("in the forceset list not found in model formula. Variable not forced."), sep = " "))
					}
				}
			}
		}
		####################
		# Logical controls #
		####################	
		if(time >= range) {
			stop("The range of simulation must be longer than shock time.")
		}
	}
	##################
	# Final warnings #
	##################	
	# Warn the opposite: in formula but not in estimation
	for(i in 2:length(all.vars(formula))) { # for all x variables
		if(!(all.vars(formula)[i] %in% diffs)) { #if its not in diffs
			if(!(all.vars(formula)[i] %in% names(lagdiffs))) {
				if(!(all.vars(formula)[i] %in% levels)) {
					if(!(all.vars(formula)[i] %in% names(lags))) {
						warning(paste(paste("Variable"), paste(all.vars(formula)[i]), paste("in the formula list not found lags, lagdiffs, levels, or differences. Variable ignored."), sep = " "))
					}
				}
			}
		}
	}
	###############
	# Create data #
	###############		
	# For X/Y: initialize list of DVs/IVs
	dv <- dvnamelist <- NULL
	ldvs <- lnumdvs <- ldvnamelist <- ldvset <- NULL # Initialize LDVs and number of lags
	lddvs <- ldnumdvs <- lddvnamelist <- lddvset <- NULL # Initialize LDDVs and number of lags
	dsiv <- lnumdsiv <- dsivnamelist <- dsivset <- NULL # Initialize shocked differenced IVs
	livs <- lnumivs <- livsnamelist <- livsset <- NULL # Ititialize lag IVs and number of lags
	lsiv <- lnumsiv <- lsivnamelist <- lsivset <- NULL # Initialize shocked lag IVs and number
	divs <- lnumdivs <- divnamelist <- divset <- NULL # Initialize differenced IVs
	ldsiv <- lnumldsiv <- ldsivnamelist <- ldsivset <- NULL # Initialize lag difference shocked IVs
	ldivs <- lnumldivs <- ldivsnamelist <- ldivsset <- NULL # Initialize lag differenced IVs
	siv <- sivnamelist <- sivset <- NULL # Initialize levels shocked IVs
	ivs <- ivsnamelist <- ivsset <- NULL # Initilialzie levels IVs
	trendvar <- nocons <- trendset <- NULL
	ec.message <- constant.message <- trend.message <- NULL
	if(ec == "TRUE") { 	# For Y: if error correction, we need differenced y
		# formula[[2]] is y
		ec.message <- "Error correction (EC) specified; dependent variable to be run in differences."
		dv <- as.data.frame(dshift(as.matrix(data[[as.character(formula[[2]])]]))) # place d.y in dataframe
		colnames(dv) <- dvnamelist <- paste("d", as.character(formula[[2]]), sep = ".") # name of Y
		assign(paste(dvnamelist), as.matrix(dv))
	} 
	else { 	# If just LDV model
		ec.message <- "Dependent variable to be run in levels."
		dv <- as.data.frame(as.matrix(data[as.character(formula[[2]])])) # place y in dataframe
		colnames(dv) <- dvnamelist <- as.character(formula[[2]])
		assign(paste(dvnamelist), as.matrix(dv))
	}
	# Lags
	if(length(lags)) {
		if(!(as.character(formula[[2]]) %in% names(lags))) { # If LDV not in lag
			if(noLDV == TRUE) {
				warning(paste("ARDL model with no LDV specified. Only finite dynamics are possible through lags of X (which might affect simulations). Are you sure you want this?"))
				### ADD exception here, maybe to deal with simulation
				lnumdvs <- 0	
			}
			else {
				lnumdvs <- 1 # we'll default to one lag
				warning("Lagged dependent variable added to model formula.")
				ldvs <- cbind(ldvs, lshift(as.matrix(data[as.character(formula[[2]])]), l = 1))
				v.name <- paste("l", 1, as.character(formula[[2]]), sep = ".")
				ldvnamelist <- c(ldvnamelist, v.name)
				assign(paste(v.name), lshift(as.matrix(data[as.character(formula[[2]])]), l = 1))
				if(constant == TRUE) { # If there is an intercept, use the mean
					ldvset <- c(ldvset, mean(as.matrix(data[as.character(formula[[2]])]), na.rm = T))
				} 
				else { # If not, set to 0
					ldvset <- c(ldvset, 0)
				}				
			}
		} # For all other lags
		for(i in 1:length(lags)) { # loop thru list
			if(names(lags)[[i]] == as.character(formula[[2]])) { # If it's y/LDV
				lnumdvs <- lags[[i]] # Get the lag numbers (only one y/LDV)
				for(o in 1:length(lags[[i]])) { # Loop thru obs in list
					ldvs <- cbind(ldvs, lshift(as.matrix(data[[names(lags[i])]]), l = lags[[i]][o]))
					v.name <- paste("l", lags[[i]][o], names(lags[i]), sep = ".")
					ldvnamelist <- cbind(ldvnamelist, v.name)
					assign(paste(v.name), lshift(as.matrix(data[[names(lags[i])]]), l = lags[[i]][o]))	
					if(constant == TRUE) { # If there is an intercept, use the mean
						ldvset <- c(ldvset, mean(as.matrix(data[as.character(formula[[2]])]), na.rm = T))
					} 
					else { # If not, set to 0
						ldvset <- c(ldvset, 0)
					}
				}
			} 
			else { # If it's an IV (not Y)
				if(names(lags)[[i]] %in% shockvar) { # If it's a shockvar
					lnumsiv <- lags[[i]]  # Get the lag numbers (only one shockvar)
					for(o in 1:length(lags[[i]])) { # Loop thru obs in list
						lsiv <- cbind(lsiv, lshift(as.matrix(data[[names(lags[i])]]), l = lags[[i]][o]))
						v.name <- paste("l", lags[[i]][o], names(lags[i]), sep = ".")
						lsivnamelist <- cbind(lsivnamelist, v.name)
						assign(paste(v.name), lshift(as.matrix(data[[names(lags[i])]]), l = lags[[i]][o]))
						if(names(lags)[i] %in% names(forceset)) {
							# If the variable is in forceset, grab its value
							lsivset <- c(lsivset, as.numeric(forceset[(names(forceset) == names(lags)[i])][1]))
						} 
						else {
							# If not, we use the mean
							lsivset <- c(lsivset, mean(data[[names(lags[i])]], na.rm = T))				
						}
					}
				} 
				else { # If it's not a shockvar
					lnumivs <- c(lnumivs, lags[[i]]) # Get the lag numbers
					for(o in 1:length(lags[[i]])) { # Loop thru obs in list
						livs <- cbind(livs, lshift(as.matrix(data[[names(lags[i])]]), l = lags[[i]][o]))
						v.name <- paste("l", lags[[i]][o], names(lags[i]), sep = ".")
						livsnamelist <- cbind(livsnamelist, v.name)
						assign(paste(v.name), lshift(as.matrix(data[[names(lags[i])]]), l = lags[[i]][o]))
						if(names(lags)[i] %in% names(forceset)) {
							# If the variable is in forceset, grab its value
							livsset <- c(livsset, as.numeric(forceset[(names(forceset) == names(lags)[i])][1]))
						} 
						else {
							# If not, we use the mean
							livsset <- c(livsset, mean(data[[names(lags[i])]], na.rm = T))				
						}
					}
				}
			}
		}
	} 	
	else { # Even if no lag specified, LDV required
		if(noLDV == TRUE) {		
			warning(paste(paste("ARDL model with no LDV specified. Only finite dynamics are possible through lags of X (which might affect simulations). Are you sure you want this?")))
			lnumdvs <- 0	
		}
		else {
			lnumdvs <- 1		# Default to one lag
			warning("Lagged dependent variable added to model formula.")
			ldvs <- cbind(ldvs, lshift(as.matrix(data[as.character(formula[[2]])]), l = 1))
			v.name <- paste("l", 1, as.character(formula[[2]]), sep = ".")
			ldvnamelist <- c(ldvnamelist, v.name)
			assign(paste(v.name), lshift(as.matrix(data[as.character(formula[[2]])]), l = 1))
			if(constant == TRUE) { # If there is an intercept, use the mean
				ldvset <- c(ldvset, mean(as.matrix(data[as.character(formula[[2]])]), na.rm = T))
			} 
			else { # If not, set to 0
				ldvset <- c(ldvset, 0)
			}
		}
	} 
	# Differences. Will be less complicated because of LDV issue and only first differences supported
	if (length(diffs)) {
		for(i in 1:length(diffs)){ # loop thru list
			if(diffs[i] %in% shockvar) { # If it's a shockvar
				dsiv <- cbind(dsiv, dshift(as.matrix(data[diffs[i]]))) 
				v.name <- paste("d.1", diffs[i], sep = ".") # > D.1 not supported
				dsivnamelist <- c(dsivnamelist, v.name)
				assign(paste(v.name), dshift(as.matrix(data[diffs[i]])))
				#if(names(diffs)[i] %in% names(forceset)) {
				#	# If the variable is in forceset, grab its value
				#	dsivset <- c(dsivset, forceset[[names(diffs)[i] == names(forceset)]])
				#} else {
				#	# If not, we set it to 0
				dsivset <- c(dsivset, 0)	
			}
			else { # for non-shock vars
				divs <- cbind(divs, dshift(as.matrix(data[diffs[i]])))
				v.name <- paste("d.1", diffs[i], sep = ".") # > D.1 not supported
				divnamelist <- c(divnamelist, v.name)
				assign(paste(v.name), dshift(as.matrix(data[diffs[i]])))
				#if(names(diffs)[i] %in% names(forceset)) {
				#	# If the variable is in forceset, grab its value
				#	divset <- c(divset, forceset[[names(diffs)[i] == names(forceset)]])
				#} else {
				#	# If not, we set it to 0
				divset <- c(divset, 0)		
			}
		}
	}		
	# Lagged differences
	if (length(lagdiffs)){
		for(i in 1:length(lagdiffs)) { # loop thru list
			if(names(lagdiffs)[[i]] == as.character(formula[[2]])) { # If it's y/LDV
				ldnumdvs <- lagdiffs[[i]] # Get the lag numbers (only one DV)
				for(o in 1:length(lagdiffs[[i]])) { # Loop thru obs in list
					lddvs <- cbind(lddvs, ldshift(as.matrix(data[[names(lagdiffs[i])]]), l = lagdiffs[[i]][o]))
					v.name <- paste("ld", lagdiffs[[i]][o], names(lagdiffs[i]), sep = ".")
					lddvnamelist <- cbind(lddvnamelist, v.name)
					assign(paste(v.name), ldshift(as.matrix(data[[names(lagdiffs[i])]]), l = lagdiffs[[i]][o]))
					# Set list to 0, LDV cannot be in forcelist
					lddvset <- c(lddvset, 0)
				}
			} 
			else { # If it's an IV (not form of Y)
				if(names(lagdiffs)[[i]] %in% shockvar) { # If it's a shockvar
					lnumldsiv <- lagdiffs[[i]] # Get the lag numbers (only one shockvar)
					for(o in 1:length(lagdiffs[[i]])){ # loop thru obs in list
						ldsiv <- cbind(ldsiv, ldshift(as.matrix(data[[names(lagdiffs[i])]]), l = lagdiffs[[i]][o]))
						v.name <- paste("ld", lagdiffs[[i]][o], names(lagdiffs[i]), sep = ".")
						ldsivnamelist <- c(ldsivnamelist, v.name)
						assign(paste(v.name), ldshift(as.matrix(data[[names(lagdiffs[i])]]), l = lagdiffs[[i]][o]))
						#if(names(lagdiffs)[i] %in% names(forceset)) {
						#	# If the variable is in forceset, grab its value
						#	ldsivset <- c(ldsivset, forceset[[names(lagdiffs)[i] == names(forceset)]])
						#} else {
						#	# If not, we set it to 0
						ldsivset <- c(ldsivset, 0)	
					}							
				} 
				else { # non shock-vars
					lnumldivs <- c(lnumldivs, lagdiffs[[i]])
					for(o in 1:length(lagdiffs[[i]])){ # loop thru obs in list
						ldivs <- cbind(ldivs, ldshift(as.matrix(data[[names(lagdiffs[i])]]), l = lagdiffs[[i]][o]))
						v.name <- paste("ld", lagdiffs[[i]][o], names(lagdiffs[i]), sep = ".")
						ldivsnamelist <- c(ldivsnamelist, v.name)
						assign(paste(v.name), ldshift(as.matrix(data[[names(lagdiffs[i])]]), l = lagdiffs[[i]][o]))
						#if(names(lagdiffs)[i] %in% names(forceset)) {
						#	# If the variable is in forceset, grab its value
						#	ldivsset <- c(ldivsset, forceset[[names(lagdiffs)[i] == names(forceset)]])
						#} else {
						#	# If not, we set it to 0
						ldivsset <- c(ldivsset, 0)		
					}
				}
			}
		}
	}
	# Levels
	if (length(levels)){
		for(i in 1:length(levels)){
			if(levels[i] %in% shockvar) { # If it's a shockvar
				siv <- cbind(siv, as.matrix(data[levels[[i]]])) 
				v.name <- levels[i]
				sivnamelist <- c(sivnamelist, v.name)
				assign(paste(v.name), as.matrix(data[levels[[i]]]))
				if(levels[i] %in% names(forceset)) {
					sivset <- c(sivset, as.numeric(forceset[levels[i]]))
				} 
				else {
					# If not, we use the mean
					sivset <- c(sivset, mean(as.matrix(data[levels[[i]]]), na.rm = T))
				}
			} 
			else {
				ivs <- cbind(ivs, as.matrix(data[levels[[i]]]))
				v.name <- levels[i]
				ivsnamelist <- c(ivsnamelist, v.name)
				assign(paste(v.name), as.matrix(data[levels[[i]]]))
				if(levels[i] %in% names(forceset)) {
					ivsset <- c(ivsset, as.numeric(forceset[levels[i]]))
				} 
				else {
					# If not, we use the mean
					ivsset <- c(ivsset, mean(as.matrix(data[levels[[i]]]), na.rm = T))
				}
			}
		}
	}	
	# Trend/constant
	if(trend == TRUE) {
		trendvar <- seq(1, length(dv[,1]), 1)
		trendset <- -burnin
		trend.message <- "Deterministic linear trend added to model formula."
	} 
	if(constant == FALSE) {
		constant.message <- "Constant suppressed from model formula."
	}	
	
	# Independent variables, their names, and the setlist of values for simulations we created
	#  Order is important here!
	IVs <- as.data.frame(cbind(ldvs, lddvs, dsiv, lsiv, ldsiv, siv, ivs, divs, livs, ldivs))
	colnames(IVs) <- c(ldvnamelist, lddvnamelist, dsivnamelist, lsivnamelist, ldsivnamelist, 
		sivnamelist, ivsnamelist, divnamelist, livsnamelist, ldivsnamelist)
	set <- setlist <- c(ldvset, lddvset, dsivset, lsivset, ldsivset, sivset, ivsset, 
		divset, livsset, ldivsset, trendset)

	# If there's a constant, add this back to the setlist, too
	if(constant == TRUE) {
		set <- setlist <- c(1, setlist)
	}
	##################
	# Estimate model #
	##################	
	if(trend == TRUE) {
		IVs <- data.frame(IVs, trendvar)
		colnames(IVs) <- c(ldvnamelist, lddvnamelist, dsivnamelist, lsivnamelist, ldsivnamelist, 
			sivnamelist, ivsnamelist, divnamelist, livsnamelist, ldivsnamelist, "trendvar")
		if(constant == FALSE) {
			res <- lm(as.formula(paste(paste(dvnamelist), "~", paste(colnames(IVs), collapse = "+"), "- 1", collapse = " ")))
		} 
		else {
			res <- lm(as.formula(paste(paste(dvnamelist), "~", paste(colnames(IVs), collapse = "+"), collapse = " ")))
		}
	} 
	else {
		if(constant == FALSE) {
			res <- lm(as.formula(paste(paste(dvnamelist), "~", paste(colnames(IVs), collapse = "+"), "- 1", collapse = " ")))
		} 
		else {
			res <- lm(as.formula(paste(paste(dvnamelist), "~", paste(colnames(IVs), collapse = "+"), collapse = " ")))
		}
	}
	if(modelout == TRUE) {
		print(summary(res))
	}
	
	# Add EC status to res object to help with pssbounds
	# Add y to help with dynardl.auto.correlated
	res$y <- dv
	res$y.name <- as.character(formula[[2]])
	res$EC <- ec
	res$simulate <- simulate
	res$ymean <- mean(as.matrix(data[[as.character(formula[[2]])]]), na.rm = T)	
	res$ldv <- ifelse(noLDV == TRUE, FALSE, TRUE)
	
	if(ec == TRUE) { # Include differenced (which should be 0) if it's EC
		res$ymean.diff <- mean(dshift(as.matrix(data[[as.character(formula[[2]])]])), na.rm = T)
	}
	
	# Warning about lags in EC model
	if(ec == TRUE) {
		if(max(lnumdvs, lnumsiv, lnumivs) > 1) {
			warning("Multiple lags included in implied co-integrated relationship: are you sure you want this?")
		}
	}
	
	# Print any messages
	print(ec.message)
	if(!(identical(NULL, trend.message))) {print(trend.message)}
	if(!(identical(NULL, constant.message))) {print(constant.message)}
	if(!(identical(NULL, shock.message))) {print(shock.message)}

	# Quantities of interest
	B <- coef(res) # get betas
	V <- vcov(res) # get vcov matrix
	sigma2 <- sigma(res)^2
	dfsig <- res$df.residual
	len <- res$rank
	
	###############
	# Simulations #
	###############		
	if(simulate == TRUE) {	
		# Set parameters
		sigl <- ((100-sig)/2)/100
		sigu <- (100-((100-sig)/2))/100
		brange <- range + burnin
		btime <- time + burnin
	
		# Vectors for predicted values and significance
		meanpv <- meandpv <- rep(NA, brange)
		# For significance: if the user sets a value of either 75, 90, or 95, we proceed as normal
		if(sig %in% c(75, 90, 95)) {
			d_PV_pctile <- PV_pctile <- matrix(rep(NA, brange*6), ncol = 6)
			colnames(PV_pctile) <- c("ll95", "ll90", "ll75", "ul75", "ul90", "ul95")
			colnames(d_PV_pctile) <- c("d.ll95", "d.ll90", "d.ll75", "d.ul75", "d.ul90", "d.ul95")
		}
		# If not, we add two columns on the outside to be custom user values
		else {
			d_PV_pctile <- PV_pctile <- matrix(rep(NA, brange*8), ncol = 8)
			colnames(PV_pctile) <- c("ll95", "ll90", "ll75", "ul75", "ul90", "ul95", "ll", "ul")
			colnames(d_PV_pctile) <- c("d.ll95", "d.ll90", "d.ll75", "d.ul75", "d.ul90", "d.ul95", "d.ll", "d.ul")
		}
		
		if(fullsims == TRUE) {
			PV_all_sims <- matrix(rep(NA, brange*sims), ncol = brange)
		}
		
		##################################
		# Values at t = 1 (first burnin) #
		##################################
		# Draw from a multivariate normal, place in PB
		PB <- mvrnorm(n = sims, mu = B, Sigma = V)
		# Replicating Stata: Sigma2 is going to be sigma2*the degrees of freedom divided by random chisq draws
		# This has length of the number of sims 
		Sigma2 <- (sigma2*dfsig)/(rchisq(sims, dfsig))
		
		# now multiply [sims x k] posterior beta matrix by [k x 1] set matrix
		PV <- PB%*%set
		
		# Expected values get error
		if(expectedval == TRUE) {
			for(i in 1:sims) { # For each PB, of which there are sims
				PV.error <- PV[i] + rnorm(1000, 0, sqrt(Sigma2[i])) # Replace it with avg of 1000 errors with a variance of the ith Sigma
				PV[i] <- mean(PV.error)
			}
		} 
		else {
			PV <- PV + rnorm(sims, 0, sqrt(Sigma2))
		}
		# Quantities depend on user values
		if(sig %in% c(75, 90, 95)) {
			# Upper/lower bounds depend on specifications
			if(ec == TRUE) { # If ECM, add percentiles to mean
				if(constant == FALSE) { # Set this because it changes where LDV is
					PV_pctile[1,] <- quantile(PV, c(0.025, 0.05, 0.125, 0.875, 0.95, 0.975)) + mean(as.matrix(data[[as.character(formula[[2]])]]), na.rm = T)
				}
				else {
					PV_pctile[1,] <- quantile(PV, c(0.025, 0.05, 0.125, 0.875, 0.95, 0.975)) + mean(as.matrix(data[[as.character(formula[[2]])]]), na.rm = T)
				}
				d_PV_pctile[1,] <- quantile(PV, c(0.025, 0.05, 0.125, 0.875, 0.95, 0.975)) # Remember, the predictions are the differences since the model is delta_y
			}
			else {  # For LDV
				PV_pctile[1,] <- quantile(PV, c(0.025, 0.05, 0.125, 0.875, 0.95, 0.975))
				d_PV_pctile[1,] <- rep(NA, 6) # the first set of differences will be empty (no period to difference)
			}
		} 
		else { # If the user needs their own significance value too
			if(ec == TRUE) { # If ECM, add percentiles to old predicted values
				if(constant == FALSE) { # Set this because it changes where LDV is
					PV_pctile[1,] <- quantile(PV, c(0.025, 0.05, 0.125, 0.875, 0.95, 0.975, sigl, sigu)) + mean(as.matrix(data[[as.character(formula[[2]])]]), na.rm = T)	
				} 
				else {
					PV_pctile[1,] <- quantile(PV, c(0.025, 0.05, 0.125, 0.875, 0.95, 0.975, sigl, sigu)) + mean(as.matrix(data[[as.character(formula[[2]])]]), na.rm = T)	
				}
				d_PV_pctile[1,] <- quantile(PV, c(0.025, 0.05, 0.125, 0.875, 0.95, 0.975, sigl, sigu)) # Remember, the predictions are the differences since the model is delta_y
			} 
			else { # For LDV
				PV_pctile[1,] <- quantile(PV, c(0.025, 0.05, 0.125, 0.875, 0.95, 0.975, sigl, sigu))
				d_PV_pctile[1,] <- rep(NA, 8) # the first set of differences will be empty (no period to difference)
			}
		}
		
		if(fullsims == TRUE) {
			if(ec == TRUE) { # Store them as predicted levels of Y
				if(constant == FALSE) { # Set this because it changes where LDV is
					PV_all_sims[,1] <- mean(PV) + mean(as.matrix(data[[as.character(formula[[2]])]]), na.rm = T)
				}
				else {
					PV_all_sims[,1] <- mean(PV) + mean(as.matrix(data[[as.character(formula[[2]])]]), na.rm = T)
				}
			} 
			else {
				PV_all_sims[,1] <- PV
			}
		}
				
		# If ECM, predicted change is added on to sample mean of Y: first/second element of forcesetlist (first is constant)
		if(ec == TRUE) {
			if(constant == FALSE) { # Set this because it changes where LDV is
				if(qoi == "mean") { # If we're summarizing with the mean
					meanpv[1] <- mean(PV) + mean(as.matrix(data[[as.character(formula[[2]])]]), na.rm = T) # For levels, we add the prediction onto the first mean value
					meandpv[1] <- mean(PV) # For differences, the prediction is the difference: the model is in delta_y
				}
				else {
					meanpv[1] <- median(PV) + mean(as.matrix(data[[as.character(formula[[2]])]]), na.rm = T) 
					meandpv[1] <- median(PV)					
				}
			} 
			else {
				if(qoi == "mean") { # If we're summarizing with the mean
					meanpv[1] <- mean(PV) + mean(as.matrix(data[[as.character(formula[[2]])]]), na.rm = T) 
					meandpv[1] <- mean(PV)						
				}
				else {
					meanpv[1] <- median(PV) + mean(as.matrix(data[[as.character(formula[[2]])]]), na.rm = T) 
					meandpv[1] <- median(PV)				
				}
			}
		} 
		else { 	# If just LDV model, it's just the new predicted level
			if(qoi == "mean") { # If we're summarizing with the mean
				meanpv[1] <- mean(PV)
				meandpv[1] <- NA # The first one will be empty: nothing to difference				
			}
			else {
				meanpv[1] <- median(PV)
				meandpv[1] <- NA # The first one will be empty: nothing to difference						
			}
		}			
		###########################
		# Values at t = 2 - range #
		###########################		
		# Progress bar to monitor dynardl
		
		print("dynardl estimating ...") 
		flush.console()
		pb <- txtProgressBar(min = 0, max = brange, style = 3)
		
		for(p in 2:brange) {
			Sys.sleep(0.1)
			PV_lag <- PV # keep the last set of PVs for differencing at the end
			row <- 1			# To work through setlist
			if(constant == TRUE) {
				row <- row + 1	# Increment past constant
			}
			
			if(max(lnumdvs) > 0) { # If LDV is in model (which it should be, if not forced not)
				##### LDVs (MUST HAVE) (these don't depend on shocktime)
				for(lag in 1:length(lnumdvs)) {		# For each lag in the LDVs 
					w <- p - lnumdvs[lag]			# w is the time period minus the appropriate lag of the LDV
					if(w > 0) {	# If LDV exists, meaning we've gone forward enough in time
						set[row] <- meanpv[w] 	# For that lag of LDV, give it the corresponding PV
					}	# Else: keep the set where it is, meaning don't replace the lag
					row <- row + 1					# Move down the setlist (i.e. to the next lag of LDV)
				}				
			}		
			##### LDDVs (OPTIONAL) (these don't depend on shocktime)			
			if(length(ldnumdvs)) {			
				for(lag in 1:length(ldnumdvs)) {	# For each lag in the LDDVs					
					w <- p - ldnumdvs[lag]			# w is the time period minus the appropriate lag of the LDV
					if(w > 1) {	# If LDDV exists, meaning we've gone forward enough in time
						wm1 <- w - 1
						set[row] <- meanpv[w] - meanpv[wm1]	# For that LDDV, give it the corresponding PV
					}	# Else: keep the set where it is, meaning don't replace the lag
					row <- row + 1					# Move down the setlist (i.e. to the next lag of LDV)
				}	
			}
			##### DSIV, LSIV, LDSIV, SIV (these depend on shocktime)
			# If it's the shocktime
			if(p == btime) { # Handles dsiv, lsiv, ldsiv, siv
				if(length(dsivnamelist)) {
					# Shock differenced var
					for(var in 1:length(dsivnamelist)) {
						set[row] <- shockval
						row <- row + 1
					}
				}
				if(length(lsivnamelist)) { # Nothing happens to LAGS in the shock period
					for(var in 1:length(lsivnamelist)) { # One var for each lag
						row <- row + 1
					}
				}
				if(length(ldsivnamelist)) { # Nothing happens to LAGGED DIFFS in the shock period
					for(var in 1:length(ldsivnamelist)) { # One var for each lag d
						row <- row + 1
					}
				}
				if(length(sivnamelist)) { # Whatever the shocked variable was in levels + shockval (only one shockvar)
					for(var in 1:length(sivnamelist)) {
						set[row] <- sivset[1] + shockval
						row <- row + 1 
					}
				}
			} 
			else { # If it's after the shocktime
				if(p > btime) { # Handles dsiv, lsiv, ldsiv, siv
					if(length(dsivnamelist)) {
						# Shock differenced var back to 0
						for(var in 1:length(dsivnamelist)) {
							set[row] <- 0
							row <- row + 1
						}
					}
					if(length(lsivnamelist)) {
						for(lag in 1:length(lnumsiv)) { # We'll use the lag list this time instead of variable names
							w <- p - btime
							if(w == lnumsiv[lag]) { # If the lag is now the ``time'' (i.e. the lag + the shocktime)
								set[row] <- lsivset[1] + shockval # Give it the mean (or forceset) plus the shock
							} # else { # Else, keep the set where it is, meaning don't replace the lag
							#set[row] <- lsivset[1]
							#}
							row <- row + 1
						}
					}
					if(length(ldsivnamelist)) {
						for(lagd in 1:length(lnumldsiv)) {
							w <- p - btime
							if(w == lnumldsiv[lagd]) { # If the lag is ``now'', the lagd variable gets the shockval
								set[row] <- shockval
							} else { # Or 0
								set[row] <- 0
							}
							row <- row + 1
						}
					}
					if(length(sivnamelist)) { # Just increment past: no return to value
						for(var in 1:length(sivnamelist)) {
							# set[row] <- sivset[1]  	# NOT returning to original value
							row <- row + 1
						}
					}		
				}
			}
			##### TREND (OPTIONAL) (doesn't depend on shocktime)
			if(trend == TRUE) {
				# It's the last coefficient in the set
				set[length(B)] <- p - burnin # trend = sim time - burnin time		
			}
			
			
			# Still in the p:brange loop. Just FYI
			# Create next predicted values. ``set'' has now been replaced by all of the new values
			PV <- PB%*%set
			# Same as before. 
			if(expectedval == TRUE) {
				for(i in 1:sims) { # For each PB, of which there are sims
					PV.error <- PV[i] + rnorm(1000, 0, sqrt(Sigma2[i])) # Average of 1000 errors with a variance of the ith Sigma
					PV[i] <- mean(PV.error)
				}
			} else {
				PV <- PV + rnorm(sims, 0, sqrt(Sigma2))
			}
			
			# Quantities depend on user values
			if(sig %in% c(75, 90, 95)) {
				# Upper/lower bounds depend on specifications
				if(ec == TRUE) { # If ECM, add percentiles to old predicted values
					if(noLDV == TRUE) { # If no LDV, then hard reference the mean, since no history through LDV
						if(constant == FALSE) { # Set this because it changes where LDV is
							PV_pctile[p,] <- quantile(PV, c(0.025, 0.05, 0.125, 0.875, 0.95, 0.975)) + mean(as.matrix(data[[as.character(formula[[2]])]]), na.rm = T)
						}
						else {
							PV_pctile[p,] <- quantile(PV, c(0.025, 0.05, 0.125, 0.875, 0.95, 0.975)) + mean(as.matrix(data[[as.character(formula[[2]])]]), na.rm = T)
						}
					}
					else { # If there is an LDV, we need to use set[] rather than the hard value			
						if(constant == FALSE) { # Set this because it changes where LDV is
							PV_pctile[p,] <- quantile(PV, c(0.025, 0.05, 0.125, 0.875, 0.95, 0.975)) + set[1]
						}
						else {
							PV_pctile[p,] <- quantile(PV, c(0.025, 0.05, 0.125, 0.875, 0.95, 0.975)) + set[2]
						}
					}
					# If it's an ECM, the PVs are the differences as the model is in delta_y
					d_PV_pctile[p,] <- quantile(PV, c(0.025, 0.05, 0.125, 0.875, 0.95, 0.975))				
				}
				else {
					PV_pctile[p,] <- quantile(PV, c(0.025, 0.05, 0.125, 0.875, 0.95, 0.975))	
					# differenced values don't depend on ECM: it's built in to the PVs
					d_PV <- PV - PV_lag
					d_PV_pctile[p,] <- quantile(d_PV, c(0.025, 0.05, 0.125, 0.875, 0.95, 0.975))
				}
			} 
			else { # If the user needs their own significance value too
				if(ec == TRUE) { # If ECM, add percentiles to old predicted values
					if(noLDV == TRUE) { # If no LDV, then hard reference the mean, since no history through LDV					
						if(constant == FALSE) { # Set this because it changes where LDV is
							PV_pctile[p,] <- quantile(PV, c(0.025, 0.05, 0.125, 0.875, 0.95, 0.975, sigl, sigu)) + mean(as.matrix(data[[as.character(formula[[2]])]]), na.rm = T)
						} 
						else {
							PV_pctile[p,] <- quantile(PV, c(0.025, 0.05, 0.125, 0.875, 0.95, 0.975, sigl, sigu)) + mean(as.matrix(data[[as.character(formula[[2]])]]), na.rm = T)
						}
					} 
					else {
						if(constant == FALSE) { # Set this because it changes where LDV is
							PV_pctile[p,] <- quantile(PV, c(0.025, 0.05, 0.125, 0.875, 0.95, 0.975, sigl, sigu)) + set[1]
						} 
						else {
							PV_pctile[p,] <- quantile(PV, c(0.025, 0.05, 0.125, 0.875, 0.95, 0.975, sigl, sigu)) + set[2]
						}
					}
					# If it's an ECM, the PVs are the differences as the model is in delta_y
					d_PV_pctile[p,] <- quantile(PV, c(0.025, 0.05, 0.125, 0.875, 0.95, 0.975, sigl, sigu))
				} 	
				else {
					PV_pctile[p,] <- quantile(PV, c(0.025, 0.05, 0.125, 0.875, 0.95, 0.975, sigl, sigu)) 
					# If it's an LDV, the PVs need to be obtained from the differences
					d_PV <- PV - PV_lag
					d_PV_pctile[p,] <- quantile(d_PV, c(0.025, 0.05, 0.125, 0.875, 0.95, 0.975, sigl, sigu))
				}
			}
			if(fullsims == TRUE) {
				if(ec == TRUE) {
					if(noLDV == TRUE) { # If no LDV, then hard reference the mean, since no history through LDV
						if(constant == FALSE) { # Set this because it changes where LDV is
							PV_all_sims[,p] <- PV + mean(as.matrix(data[[as.character(formula[[2]])]]), na.rm = T)
						}
						else {
							PV_all_sims[,p] <- PV + mean(as.matrix(data[[as.character(formula[[2]])]]), na.rm = T) 
						}
					} 
					else {
						if(constant == FALSE) { # Set this because it changes where LDV is
							PV_all_sims[,p] <- PV + set[1]
						}
						else {
							PV_all_sims[,p] <- PV + set[2]
						}
					} 
				}
				else {
					PV_all_sims[,p] <- PV
				}
			}
			
			# Lastly, get the predicted values
			if(ec == TRUE) {
				if(constant == FALSE) {
					if(qoi == "mean") { # If we're summarizing with the mean
						if(noLDV == TRUE) { # If no LDV, then hard reference the mean, since no history through LDV
							meanpv[p] <- mean(PV) + mean(as.matrix(data[[as.character(formula[[2]])]]), na.rm = T)
						} 
						else {
							meanpv[p] <- mean(PV) + set[1]							
						}
						meandpv[p] <- mean(PV) # For differences, the prediction is the difference: the model is in delta_y
					}
					else { # if it's the median
						if(noLDV == TRUE) { # If no LDV, then hard reference the mean, since no history through LDV
							meanpv[p] <- median(PV) + mean(as.matrix(data[[as.character(formula[[2]])]]), na.rm = T)
						}
						else {
							meanpv[p] <- median(PV) + set[1]							
						}
						meandpv[p] <- median(PV) # For differences, the prediction is the difference: the model is in delta_y
					}
				} 
				else {
					if(qoi == "mean") { # If we're summarizing with the mean
						if(noLDV == TRUE) { # If no LDV, then hard reference the mean, since no history through LDV
							meanpv[p] <- mean(PV) + mean(as.matrix(data[[as.character(formula[[2]])]]), na.rm = T)
						}
						else {
							meanpv[p] <- mean(PV) + set[2]							
						}
						meandpv[p] <- mean(PV) # For differences, the prediction is the difference: the model is in delta_y
					}
					else { # if it's the median
						if(noLDV == TRUE) { # If no LDV, then hard reference the mean, since no history through LDV
							meanpv[p] <- median(PV) + mean(as.matrix(data[[as.character(formula[[2]])]]), na.rm = T)
						}
						else { # if it's the median
							meanpv[p] <- median(PV) + set[2]							
						}
						meandpv[p] <- median(PV)	 # For differences, the prediction is the difference: the model is in delta_y
					}
				}
			} 	
			else { # if it's in levels
				if(qoi == "mean") { # If we're summarizing with the mean
					meanpv[p] <- mean(PV)
					meandpv[p] <- mean(d_PV)			
				}
				else {
					meanpv[p] <- median(PV)		
					meandpv[p] <- median(d_PV)								
				}
			}				
			setTxtProgressBar(pb, p)
		} # Close time loop
		close(pb)		
		
		###################################
		# Get simulation values to report #
		###################################
		# Discard the burnins
		if(sig %in% c(75, 90, 95)) {
			sims <- matrix(rep(NA, range*14), ncol = 14)
			sims[,1] <- meanpv[(burnin+1):brange]
			sims[,2:7] <- PV_pctile[(burnin+1):brange,]
			sims[,8] <- meandpv[(burnin+1):brange]
			sims[,9:14] <- d_PV_pctile[(burnin+1):brange,]
			colnames(sims) <- c("central", "ll95", "ll90", "ll75", "ul75", "ul90", "ul95", 
				"d.central",  "d.ll95", "d.ll90", "d.ll75", "d.ul75", "d.ul90", "d.ul95")
		} 
		else {
			sims <- matrix(rep(NA, range*18), ncol = 18)
			sims[,1] <- meanpv[(burnin+1):brange]
			sims[,2:9] <- PV_pctile[(burnin+1):brange,]
			sims[,10] <- meandpv[(burnin+1):brange]
			sims[,11:18] <- d_PV_pctile[(burnin+1):brange,]
			colnames(sims) <- c("central", "ll95", "ll90", "ll75", "ul75", "ul90", "ul95", paste("ll", sig, sep = ""), paste("ul", sig, sep = ""),
				"d.central", "d.ll95", "d.ll90", "d.ll75", "d.ul75", "d.ul90", "d.ul95", paste("d.ll", sig, sep = ""), paste("d.ul", sig, sep = ""))
		}
		sim.time <- seq(1, length(sims[,1]))
		temp.names <- colnames(sims)
		sims <- cbind(sim.time, sims)
		colnames(sims) <- c("time", temp.names)
		z <- data.frame(sims)
		
		if(fullsims == TRUE) {
			all.sims <- data.frame(PV_all_sims[,(burnin+1):brange])
			colnames(all.sims) <- paste("time", seq(1, ncol(all.sims), 1), sep = "")
			all.sims$central <- paste(qoi)
		}		
		
		#########################
		# Establish data output #
		#########################
		if(fullsims == TRUE) {
				out <- list(z, res, all.sims)
				names(out) <- c("simulation", "model", "rawsims")
				out$simulation$shocktime <- time				
		}
		else {
			out <- list(z, res)
			names(out) <- c("simulation", "model")
			out$simulation$shocktime <- time
		}
	} # This closes the simulation brackets
	else { # If simulation is false
		out <- list(res)
		names(out) <- c("model")
	}
	class(out) <- "dynardl"  # Apply class dynardl so summary method is called (and for pssbounds)
	out	# Send data out	
}

#########################################
# ------(5) area.simulation.plot -------#
#########################################
#' Create an area plot of a simulated response in a dynardl model
#' @param x a dynardl model with a simulation to be plotted
#' @param response whether the plot of the response should be shown in levels of the dependent variable (\code{levels}) 
#' or in changes from the mean of the dependent variable (\code{mean.changes}). The default is \code{levels}
#' @param bw should the colors be in black and white (for publication)? The default is \code{FALSE}
#' @return an area plot
#'
#' @name dynamac-deprecated
#' @seealso \code{link{dynamac-deprecated}}
#' @keywords internal
NULL

#' @rdname dynamac-deprecated
#' @section \code{area.simulation.plot}:
#' For \code{area.simulation.plot}, use \code{\link{dynardl.simulation.plot}}.
#' 
#' @importFrom graphics lines plot points polygon segments
#' @author Soren Jordan and Andrew Q. Philips
#' @keywords utilities
#' @export

area.simulation.plot <- function(x, response = "levels", bw = FALSE) {
	.Deprecated("dynardl.simulation.plot")
	if(x$model$simulate == FALSE) {
		stop("dynardl object does not include simulation to plot.")
	}
	if(!(response %in% c("levels", "mean.changes"))) {
		stop("Response must either be plotted in levels ('levels') or in changes from mean of DV ('mean.changes').")
	}
	z <- x$simulation
	if(response == "levels") { # If we're just plotting levels of Y
		plot(z$time, z$ll95, type = "n", ylim = c(min(z$ll95), max(z$ul95)), 
			ylab = "Y Value", xlab = "Time")
		if(bw == FALSE) {
			# 95
			polygon(c(z$time, rev(z$time)), c(z$ul95, rev(z$ll95)), col = "skyblue1", border = NA)
			# 90
			polygon(c(z$time, rev(z$time)), c(z$ul90, rev(z$ll90)), col = "skyblue3", border = NA)
			# 75
			polygon(c(z$time, rev(z$time)), c(z$ul75, rev(z$ll75)), col = "grey30", border = NA)
			# Actual response
			lines(z$time, z$central, lty = 2, lwd = 3)
		} else {
			# 95
			polygon(c(z$time, rev(z$time)), c(z$ul95, rev(z$ll95)), col = "grey70", border = NA)
			# 90
			polygon(c(z$time, rev(z$time)), c(z$ul90, rev(z$ll90)), col = "grey50", border = NA)
			# 75
			polygon(c(z$time, rev(z$time)), c(z$ul75, rev(z$ll75)), col = "grey30", border = NA)
			# Actual response
			lines(z$time, z$central, lty = 2, lwd = 3)	
		}	
	} else { # If it's changes from the mean, changes values, same code
		z$ll95 <- z$ll95 - x$model$ymean
		z$ul95 <- z$ul95 - x$model$ymean
		z$ll90 <- z$ll90 - x$model$ymean
		z$ul90 <- z$ul90 - x$model$ymean
		z$ll75 <- z$ll75 - x$model$ymean
		z$ul75 <- z$ul75 - x$model$ymean
		z$central <- z$central - x$model$ymean
		plot(z$time, z$ll95, type = "n", ylim = c(min(z$ll95), max(z$ul95)), 
			ylab = "Changes from Y Mean Value", xlab = "Time")
		if(bw == FALSE) {
			# 95
			polygon(c(z$time, rev(z$time)), c(z$ul95, rev(z$ll95)), col = "skyblue1", border = NA)
			# 90
			polygon(c(z$time, rev(z$time)), c(z$ul90, rev(z$ll90)), col = "skyblue3", border = NA)
			# 75
			polygon(c(z$time, rev(z$time)), c(z$ul75, rev(z$ll75)), col = "grey30", border = NA)
			# Actual response
			lines(z$time, z$central, lty = 2, lwd = 3)
		} else {
			# 95
			polygon(c(z$time, rev(z$time)), c(z$ul95, rev(z$ll95)), col = "grey70", border = NA)
			# 90
			polygon(c(z$time, rev(z$time)), c(z$ul90, rev(z$ll90)), col = "grey50", border = NA)
			# 75
			polygon(c(z$time, rev(z$time)), c(z$ul75, rev(z$ll75)), col = "grey30", border = NA)
			# Actual response
			lines(z$time, z$central, lty = 2, lwd = 3)	
		}	
	}
}

#########################################
# ------(6) spike.simulation.plot ------#
#########################################
#' Create a spike plot of a simulated response in a dynardl model
#' @param x a dynardl model with a simulation to be plotted
#' @param response whether the plot of the response should be shown in levels of the dependent variable (\code{levels}) 
#' or in changes from the mean of the dependent variable (\code{mean.changes}). The default is \code{levels}
#' @param bw should the colors be in black and white (for publication)? The default is \code{FALSE}
#' @return a spike plot
#'
#' @name dynamac-deprecated
#' @seealso \code{link{dynamac-deprecated}}
#' @keywords internal
NULL

#' @rdname dynamac-deprecated
#' @section \code{spike.simulation.plot}:
#' For \code{spike.simulation.plot}, use \code{\link{dynardl.simulation.plot}}.
#'
#' @importFrom graphics lines plot points polygon segments
#' @author Soren Jordan and Andrew Q. Philips
#' @keywords utilities
#' @export

spike.simulation.plot <- function(x, response = "levels", bw = FALSE) {
	.Deprecated("dynardl.simulation.plot")
	if(x$model$simulate == FALSE) {
		stop("dynardl object does not include simulation to plot.")
	}
	if(!(response %in% c("levels", "mean.changes"))) {
		stop("Response must either be plotted in levels ('levels') or in changes from mean of DV ('mean.changes').")
	}
	z <- x$simulation
	if(response == "levels") { # If we're just plotting levels of Y
		plot(z$time, z$ll95, type = "n", ylim = c(min(z$ll95), max(z$ul95)), 
			ylab = "Y Value", xlab = "Time")	
		if(bw == FALSE) {
			for(i in 1:length(z$time)) { # 95 percent sig
				segments(z$time[i], z$ll95[i], z$time[i], z$ul95[i], lwd = 1, col = "skyblue1")
			}
			for(i in 1:length(z$time)) { # 90 percent sig
				segments(z$time[i], z$ll90[i], z$time[i], z$ul90[i], lwd = 3, col = "skyblue3")
			}
			for(i in 1:length(z$time)) { # 75 percent sig
				segments(z$time[i], z$ll75[i], z$time[i], z$ul75[i], lwd = 5, col = "grey30")
			}
			# Actual response
			points(z$time, z$central, lwd = 4)
		} else {
			for(i in 1:length(z$time)) { # 95 percent sig
				segments(z$time[i], z$ll95[i], z$time[i], z$ul95[i], lwd = 1, col = "grey70")
			}
			for(i in 1:length(z$time)) { # 90 percent sig
				segments(z$time[i], z$ll90[i], z$time[i], z$ul90[i], lwd = 3, col = "grey50")
			}
			for(i in 1:length(z$time)) { # 75 percent sig
				segments(z$time[i], z$ll75[i], z$time[i], z$ul75[i], lwd = 5, col = "grey30")
			}
			# Actual response
			points(z$time, z$central, lwd = 4)	
		}
	} else {	 # If it's changes from the mean, changes values, same code
		z$ll95 <- z$ll95 - x$model$ymean
		z$ul95 <- z$ul95 - x$model$ymean
		z$ll90 <- z$ll90 - x$model$ymean
		z$ul90 <- z$ul90 - x$model$ymean
		z$ll75 <- z$ll75 - x$model$ymean
		z$ul75 <- z$ul75 - x$model$ymean
		z$central <- z$central - x$model$ymean
		plot(z$time, z$ll95, type = "n", ylim = c(min(z$ll95), max(z$ul95)), 
			ylab = "Changes from Y Mean Value", xlab = "Time")	
		if(bw == FALSE) {
			for(i in 1:length(z$time)) { # 95 percent sig
				segments(z$time[i], z$ll95[i], z$time[i], z$ul95[i], lwd = 1, col = "skyblue1")
			}
			for(i in 1:length(z$time)) { # 90 percent sig
				segments(z$time[i], z$ll90[i], z$time[i], z$ul90[i], lwd = 3, col = "skyblue3")
			}
			for(i in 1:length(z$time)) { # 75 percent sig
				segments(z$time[i], z$ll75[i], z$time[i], z$ul75[i], lwd = 5, col = "grey30")
			}
			# Actual response
			points(z$time, z$central, lwd = 4)
		} else {
			for(i in 1:length(z$time)) { # 95 percent sig
				segments(z$time[i], z$ll95[i], z$time[i], z$ul95[i], lwd = 1, col = "grey70")
			}
			for(i in 1:length(z$time)) { # 90 percent sig
				segments(z$time[i], z$ll90[i], z$time[i], z$ul90[i], lwd = 3, col = "grey50")
			}
			for(i in 1:length(z$time)) { # 75 percent sig
				segments(z$time[i], z$ll75[i], z$time[i], z$ul75[i], lwd = 5, col = "grey30")
			}
			# Actual response
			points(z$time, z$central, lwd = 4)	
		}
	}
}	

##########################################
# ------------(7) pssbounds -------------#
##########################################
#' Perform Pesaran, Shin, and Smith (2001) cointegration test
#' @param data an optional \code{\link{dynardl}} model. This option is highly recommended. Users are welcome to supply their own case, k regressors, t-statistic, F-statistic, and observations, but it is easier to have the model determine these quantities. If a \code{\link{dynardl}} model is supplied, user-supplied arguments are ignored
#' @param obs number of observations
#' @param fstat F-statistic of the joint test that variables in first lags are equal to zero: the specific restriction tested 
#' is \code{l.y + l.1.x1 + l.1.x2 + ... + l.1.xk = 0}, except in cases II and IV (see \code{restriction} and \code{case})
#' @param tstat t-statistic of the lagged dependent variable
#' @param case The case of the test, as per Pesaran, Shin, and Smith (2001). Case I: no intercept or trend; case II: restricted intercept, no trend; case III: unrestricted intercept with no trend; case IV: unrestricted intercept and restricted trend; case V: unrestricted intercept and trend. Case III is most frequently specified
#' @param restriction if you design to test case II or IV of pssbounds, where it is assumed that the constant (case 2) or trend (case 4) are restricted in the resulting F-test, indicate that restriction = \code{TRUE}. If restriction = \code{TRUE} and there is no trend in the regression (trend = \code{FALSE} in \code{\link{dynardl}}), the F-test will include the constant in addition to the lagged dependent variable and lagged regressors in order to test for cointegration under the assumption of a restricted constant (see Pesaran, Shin and Smith [2001], case II). If restriction = \code{TRUE} and there is a trend in the regression (trend = \code{TRUE} in \code{\link{dynardl}}), the F-test will include the trend term in addition to the lagged dependent variable and lagged regressors in order to test for cointegration under the assumption of a restricted trend (see Pesaran, Shin and Smith [2001], case IV). If you are estimating the regular unrestricted ECM (this is more common), restriction = \code{FALSE}. The default is \code{FALSE}
#' @param k number of regressors appearing in levels in the estimated model, not including the lagged dependent variable
#' @param digits the number of digits to round to when showing output. The default is \code{3}
#' @param object.out if \code{TRUE}, and \code{pssbounds} is assigned to an object, the test quantities will be stored for the user's convenience
#' @details
#' pssbounds performs post-estimation cointegration testing using the bounds testing procedure from Pesaran, Shin, and Smith (2001). Since test statistics vary based on the number of \code{k} regressors, length of the series, these are required, in addition to F- and t-statistics
#' @importFrom stats coef vcov
#' @author Soren Jordan and Andrew Q. Philips
#' @keywords cointegration
#' @examples
#' # Using the ineq data from dynamac
#' # We can get all the values by hand
#' ardl.model <- dynardl(concern ~ incshare10 + urate, data = ineq, 
#'         lags = list("concern" = 1, "incshare10" = 1),
#'         diffs = c("incshare10", "urate"), 
#'         lagdiffs = list("concern" = 1),
#'         ec = TRUE, simulate = FALSE)
#' summary(ardl.model)
#' pssbounds(obs = 47, fstat = 7.01578, tstat = -3.223, case = 3, k = 1)
#' 
#' # Or just pass a dynardl model.
#' pssbounds(ardl.model)
#' @export

pssbounds <- function(data = list(), obs = NULL, fstat = NULL, tstat = NULL, case = NULL, k = NULL, restriction = FALSE, digits = 3, object.out = FALSE) {
	# First check: was it passed something in data?
	if(identical(data, list())) { # If not (if we're working from user args)
		if(identical(obs, NULL)) {
			stop("Provide either number of observations (obs = ) or data = dynardl model object.")
		} else if(identical(fstat, NULL)) {
			stop("Provide either fstat on lagged variables (fstat = ) or data = dynardl model object.")
		} else if(identical(case, NULL)) {
			stop("Provide either case of regression (case = ) or data = dynardl model object.")
		} else if(identical(k, NULL)) {
			stop("Provide either k number of lagged variables (k = ) or data = dynardl model object.")
		}
	}
	else { # If data is not empty
		# Stops for invalid models/arguments
		if(!(identical(class(data), "dynardl"))) { # Check if data are from dynardl
			stop("To execute pssbounds, supply either a dynardl model or each argument of pssbounds.")
		} else if(data$model$EC == FALSE) {
			stop("pssbounds only for error-correcting relationships. See Philips (2018).")
		} else if(data$model$ldv == FALSE) {
			stop("pssbounds only for models that include a lagged dependent variable (LDV). See Philips (2018).")
		}	
		# Only possible if a variable (implied error-correcting) is in first lags, Isolate variables in model in first lags
		# Define three old arguments from legacy code (res, constant, and trend) from when this was in dynardl
		res <- data$model
		B <- coef(res)
		V <- vcov(res)
		constant <- ifelse("(Intercept)" %in% names(coef(res)), TRUE, FALSE)
		trend <- ifelse("trendvar" %in% names(coef(res)), TRUE, FALSE)
		R.temp <- rep(NA, length(coef(res)))
		for(i in 1:length(R.temp)) {
			# If that variable is a first lag, it's in the restriction set. 0 if not
			ifelse(grepl("l.1.", names(coef(res))[i]) == TRUE, R.temp[i] <- 1, R.temp[i] <- 0)
			# Except if it's the LDV, we don't want it
			# ifelse(names(coef(res))[i] == paste(ldvnamelist), R.temp[i] <- 0, R.temp[i] <- R.temp[i])
		}
		k.temp.fun <- sum(R.temp) - 1 # The LDV is automatically added. So if it's the only variable in first lags, this will be zero
		if(k.temp.fun < 1) { # Either if no other variables in lags, or no LDV, or both
			stop("pssbounds not executed: no variables in first lags, so no cointegrating relationship implied.")
		} 
		else {
			if(trend == TRUE) {
				if(constant == FALSE) {
					case.fun <- 0
				}
				else {
					if(restriction == TRUE) {
						case.fun <- 4
						for(i in 1:length(R.temp)) {
							ifelse(grepl("trendvar", names(coef(res))[i]) == TRUE, R.temp[i] <- 1, R.temp[i] <- R.temp[i])	
						}
					}
					else {
						case.fun <- 5
					}
				}
			}
			else {
				if(constant == FALSE) {
					case.fun <- 1
				}
				else {
					if(restriction == TRUE) {
						case.fun <- 2
						for(i in 1:length(R.temp)) {
							ifelse(grepl("(Intercept)", names(coef(res))[i]) == TRUE, R.temp[i] <- 1, R.temp[i] <- R.temp[i])	
						}
					} 
					else {
						case.fun <- 3
					}
				}
			}
			k.fun <- sum(R.temp)			
			tstat.fun <- coef(summary(res))[paste(paste("l.1.", res$y.name, sep = "")),3] # t stat on LDV
			obs.fun <- length(res$residuals) # number of observations in model
			# This needs to be expanded so that each hypothesis gets its own row
			R <- matrix(rep(0, length(coef(res))*k.fun), nrow = k.fun)
			the.row <- 1
			for(i in 1:length(R.temp)){
				if(R.temp[i] == 1) { # If that model coefficient is a first lag
					R[the.row, i] <- 1	# Then that row gets a 1 where the coefficient is
					the.row <- the.row + 1
				}			
			}
			# Restriction is always that it's equal to 0
			q <- 0
			fstat.fun <- (1/k.fun)*t(R%*%B-q)%*%solve(R%*%V%*%t(R))%*%(R%*%B-q)	
			pssbounds <- data.frame(obs.fun, k.temp.fun, tstat.fun, fstat.fun, case.fun) # k statistic here EXCLUDES the LDV
			names(pssbounds) <- c("obs", "k", "tstat", "fstat", "case")
		}
	} # End creating args from dynardl object
	# Now check against user args
	if(!(identical(data, list()))) { # If there was a model provided
		if(!(identical(NULL, obs))) { # If the user passed obs in addition to model
			if(!(identical(obs, pssbounds$obs))) { # If it doens't match dynardl, warning!
				warning(paste(paste(paste("Observations supplied ("), paste(obs), paste(") different from observations calculated by dynardl ("), 
					paste(pssbounds$obs), paste("). 	Data from dynardl used."), sep = ""), "\n", 
				paste("To execute bounds test with user data, run pssbounds() and supply each argument without a dynardl model.")))		
			}
		}
		if(!(identical(NULL, fstat))) { # If the user passed fstat in addition to model
			if(!(identical(fstat, pssbounds$fstat))) {
				warning(paste(paste(paste("F-stat supplied ("), paste(fstat), paste(") different from F-stat calculated by dynardl ("), 
					paste(pssbounds$fstat), paste("). Data from dynardl used."), sep = ""), "\n",
			paste("To execute bounds test with user data, run pssbounds() and supply each argument without a dynardl model.")))	
			}
		}
		if(!(identical(NULL, tstat))) { # If the user passed tstat in addition to model
			if(!(identical(tstat, pssbounds$tstat))) {
				warning(paste(paste(paste("t-stat supplied ("), paste(fstat), paste(") different from t-stat calculated by dynardl ("), 
					paste(pssbounds$tstat), paste("). Data from dynardl used."), sep = ""), "\n",
				paste("To execute bounds test with user data, run pssbounds() and supply each argument without a dynardl model.")))	
			}
		}
		if(!(identical(NULL, case))) { # If the user passed case in addition to model
			if(!(identical(case, pssbounds$case))) {
				warning(paste(paste(paste("Case supplied ("), paste(case), paste(") different from case calculated by dynardl ("), 
					paste(pssbounds$case), paste("). Data from dynardl used."), sep = ""), "\n",
				paste("To execute bounds test with user data, run pssbounds() and supply each argument without a dynardl model.")))
			}
		}
		if(!(identical(NULL, k))) { # If the user passed k in addition to model
			if(!(identical(k, pssbounds$k))) {
				warning(paste(paste(paste("k supplied ("), paste(case), paste(") different from k calculated by dynardl ("), 
					paste(pssbounds$k), paste("). Data from dynardl used."), sep = ""), "\n",
				paste("To execute bounds test with user data, run pssbounds() and supply each argument without a dynardl model.")))	
			}
		}
		# And finally, assign them the dynardl values anyway (after warning)
		obs <- pssbounds$obs
		fstat <- pssbounds$fstat
		tstat <- pssbounds$tstat
		case <- pssbounds$case
		k <- pssbounds$k
	}
	# Find critical values
	cases <- seq(1, 5, 1)
	cases.roman <- c("I", "II", "III", "IV", "V")
	fnote <- tnote <- cnote <- NULL
	if((case %in% cases.roman) == TRUE) {
		if(case == "I") {case <- 1}
		else if (case == "II") {case <- 2}
		else if (case == "III") {case <- 3}
		else if (case == "IV") {case <- 4}
		else if (case == "V") {case <- 5}
		}
	if((case %in% cases) == FALSE) {
		stop("Case must be 1, 2, 3, 4, or 5.")
	}
	tstat <- round(tstat, digits = digits)
	fstat <- round(fstat, digits = digits)
	if (obs <= 30) {
		# Case 1: no intercept, no trend
		if (case == 1) {
			# Doesn't exist! Use asym. crit values
			fnote <- "Small-sample critical values not provided for Case I. Asymptotic critical values used."
			#		0.10        0.05	   0.010
			# 	 I(0)   I(1)  I(0)  I(1)  I(0) I(1)
			fmat <- matrix(c(3.00, 3.00, 4.20, 4.20, 7.17, 7.17, # 0
				2.44, 3.28, 3.15, 4.11, 4.81, 6.02, # 1
				2.17, 3.19, 2.72, 3.83, 3.88, 5.30, # 2
				2.01, 3.10, 2.45, 3.63, 3.42, 4.84, # 3
				1.90, 3.01, 2.26, 3.48, 3.07, 4.44, # 4
				1.81, 2.93, 2.14, 3.34, 2.82, 4.21, # 5
				1.75, 2.87, 2.04, 3.24, 2.66, 4.05, # 6
				1.70, 2.83, 1.97, 3.18, 2.54, 3.91, # 7
				1.66, 2.79, 1.91, 3.11, 2.45, 3.79, # 8
				1.63, 2.75, 1.86, 3.05, 2.34, 3.68, # 9
				1.60, 2.72, 1.82, 2.99, 2.26, 3.60), ncol = 6, byrow = T) # 10
			if (!is.null(tstat)) {
				tnote <- "Small-sample critical values not provided for Case I. Asymptotic critical values used."
				#		0.10         0.05		   0.010
				# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
				tmat <- matrix(c(-1.62, -1.62, -1.95, -1.95, -2.58, -2.58, # 0
					-1.62, -2.28, -1.95, -2.60, -2.58, -3.22, # 1
					-1.62, -2.68, -1.95, -3.02, -2.58, -3.66, # 2
					-1.62, -3.00, -1.95, -3.33, -2.58, -3.97, # 3
					-1.62, -3.26, -1.95, -3.60, -2.58, -4.23, # 4
					-1.62, -3.49, -1.95, -3.83, -2.58, -4.44, # 5
					-1.62, -3.70, -1.95, -4.04, -2.58, -4.67, # 6
					-1.62, -3.90, -1.95, -4.23, -2.58, -4.88, # 7
					-1.62, -4.09, -1.95, -4.43, -2.58, -5.07, # 8
					-1.62, -4.26, -1.95, -4.61, -2.58, -5.25, # 9
					-1.62, -4.42, -1.95, -4.76, -2.58, -5.44), ncol = 6, byrow = T) # 10
			}
		}
		# Case II: Restricted intercept and no trend
		else if (case == 2) {
			#		0.10         0.05		   0.010
			# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
			fmat <- matrix(c(4.025, 4.025, 5.070, 5.070, 7.595, 7.595, # 0
				3.303, 3.797, 4.090, 4.663, 6.027, 6.760, # 1
				2.915, 3.695, 3.538, 4.428, 5.155, 6.265, # 2
				2.676, 3.586, 3.272, 4.306, 4.614, 5.966, # 3
				2.525, 3.560, 3.058, 4.223, 4.280, 5.840, # 4
				2.407, 3.517, 2.910, 4.193, 4.134, 5.761, # 5
				2.334, 3.515, 2.794, 4.148, 3.976, 5.691, # 6
				2.277, 3.498, 2.730, 4.163, 3.864, 5.694), ncol = 6, byrow = T) # 7
			if (!is.null(tstat)) {
				tnote <- "Critical values do not currently exist for Case II."
			}
		}
		# Case III: unrestricted intercept, no trend
		else if (case == 3)	{
			#		0.10         0.05		   0.010
			# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
			fmat <- matrix(c(6.840, 6.840, 8.770, 8.770, 13.680, 13.680, # 0
				4.290, 5.080, 5.395, 6.350, 8.170, 9.285, # 1
				3.437, 4.470, 4.267, 5.473, 6.183, 7.873, # 2
				3.008, 4.150, 3.710, 5.018, 5.333, 7.063, # 3
				2.752, 3.994, 3.354, 4.774, 4.768, 6.670, # 4
				2.578, 3.858, 3.125, 4.608, 4.537, 6.370, # 5
				2.457, 3.797, 2.970, 4.499, 4.270, 6.211, # 6
				2.384, 3.728, 2.875, 4.445, 4.104, 6.151), ncol = 6, byrow = T) # 7
			if (!is.null(tstat)) {
				tnote <- "Small-sample critical values not provided for Case III. Asymptotic critical values used."
				#		0.10         0.05		   0.010
				# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
				tmat <- matrix(c(-2.57, -2.57, -2.86, -2.86, -3.43, -3.43, # 0
					-2.57, -2.91, -2.86, -3.22, -3.43, -3.82, # 1
					-2.57, -3.21, -2.86, -3.53, -3.43, -4.10, # 2
					-2.57, -3.46, -2.86, -3.78, -3.43, -4.37, # 3
					-2.57, -3.66, -2.86, -3.99, -3.43, -4.60, # 4
					-2.57, -3.86, -2.86, -4.19, -3.43, -4.79, # 5
					-2.57, -4.04, -2.86, -4.38, -3.43, -4.99, # 6
					-2.57, -4.23, -2.86, -4.57, -3.43, -5.19, # 7
					-2.57, -4.40, -2.86, -4.72, -3.43, -5.37, # 8
					-2.57, -4.56, -2.86, -4.88, -3.42, -5.54, # 9
					-2.57, -4.69, -2.86, -5.03, -3.43, -5.68), ncol = 6, byrow = T)	# 10
			}
		}
		# Case IV: Unrestricted intercept, restricted trend
		else if (case == 4)	{
			#		0.10         0.05		   0.010
			# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
			fmat <- matrix(c(5.785, 5.785, 7.040, 7.040, 10.200, 10.200, # 0
				4.427, 4.957, 5.377, 5.963, 7.593, 8.350, # 1
				3.770, 4.535, 4.535, 5.415, 6.428, 7.505, # 2
				3.378, 4.274, 4.048, 5.090, 5.666, 6.988, # 3
				3.097, 4.118, 3.715, 4.878, 5.205, 6.640, # 4
				2.907, 4.010, 3.504, 4.743, 4.850, 6.473, # 5
				2.781, 3.941, 3.326, 4.653, 4.689, 6.358, # 6
				2.681, 3.887, 3.194, 4.604, 4.490, 6.328), ncol = 6, byrow = T) # 7
			if (!is.null(tstat)) {
				tnote <- "Critical values do not currently exist for Case IV."
			}
		}
		# Case V: unrestricted intercept, unrestricted trend
		else { # case == 5
			#		0.10         0.05		   0.010
			# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
			fmat <- matrix(c(10.340, 10.340, 12.740, 12.740, 18.560, 18.560, # 0
				6.010, 6.780, 7.360, 8.265, 10.605, 11.650, # 1
				4.577, 5.600, 5.550, 6.747, 7.977, 9.413, # 2
				3.868, 4.965, 4.683, 5.980, 6.643, 8.313, # 3
				3.430, 4.624, 4.154, 5.540, 5.856, 7.578, # 4
				3.157, 4.412, 3.818, 5.253, 5.347, 7.242, # 5
				2.977, 4.260, 3.576, 5.065, 5.046, 6.930, # 6
				2.483, 4.160, 3.394, 4.939, 4.779, 6.821), ncol = 6, byrow = T)	# 7
			if (!is.null(tstat)) {
				tnote <- "Small-sample critical values not provided for Case V. Asymptotic critical values used."
				#		0.10         0.05		   0.010
				# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
				tmat <- matrix(c(-3.13, -3.13, -3.41, -3.41, -3.96, -3.97, # 0
					-3.13, -3.40, -3.41, -3.69, -3.96, -4.26, # 1
					-3.13, -3.63, -3.41, -3.95, -3.96, -4.53, # 2
					-3.13, -3.84, -3.41, -4.16, -3.96, -4.73, # 3
					-3.13, -4.04, -3.41, -4.36, -3.96, -4.96, # 4
					-3.13, -4.21, -3.41, -4.52, -3.96, -5.13, # 5
					-3.13, -4.37, -3.41, -4.69, -3.96, -5.31, # 6
					-3.13, -4.53, -3.41, -4.85, -3.96, -5.49, # 7
					-3.13, -4.68, -3.41, -5.01, -3.96, -5.65, # 8
					-3.13, -4.82, -3.41, -5.15, -3.96, -5.79, # 9
					-3.13, -4.96, -3.41, -5.29, -3.96, -5.94), ncol = 6, byrow = T) # 10
			}
		}
	}
	else if (obs <= 35) {
		# Case I: no intercept, no trend
		if (case == 1) {
			# Doesn't exist! Use asym. crit values
			fnote <- "Small-sample critical values are not provided for Case I. Asymptotic critical values used."
			#		0.10        0.05	   0.010
			# 	 I(0)   I(1)  I(0)  I(1)  I(0) I(1)
			fmat <- matrix(c(3.00, 3.00, 4.20, 4.20, 7.17, 7.17, # 0
				2.44, 3.28, 3.15, 4.11, 4.81, 6.02, # 1
				2.17, 3.19, 2.72, 3.83, 3.88, 5.30, # 2
				2.01, 3.10, 2.45, 3.63, 3.42, 4.84, # 3
				1.90, 3.01, 2.26, 3.48, 3.07, 4.44, # 4
				1.81, 2.93, 2.14, 3.34, 2.82, 4.21, # 5
				1.75, 2.87, 2.04, 3.24, 2.66, 4.05, # 6
				1.70, 2.83, 1.97, 3.18, 2.54, 3.91, # 7
				1.66, 2.79, 1.91, 3.11, 2.45, 3.79, # 8
				1.63, 2.75, 1.86, 3.05, 2.34, 3.68, # 9
				11.60, 2.72, 1.82, 2.99, 2.26, 3.60), ncol = 6, byrow = T) # 10
			if (!is.null(tstat)) {
				tnote <- "Small-sample critical values not provided for Case I. Asymptotic critical values used."
				#		0.10         0.05		   0.010
				# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
				tmat <- matrix(c(-1.62, -1.62, -1.95, -1.95, -2.58, -2.58, # 0
					-1.62, -2.28, -1.95, -2.60, -2.58, -3.22, # 1
					-1.62, -2.68, -1.95, -3.02, -2.58, -3.66, # 2
					-1.62, -3.00, -1.95, -3.33, -2.58, -3.97, # 3
					-1.62, -3.26, -1.95, -3.60, -2.58, -4.23, # 4
					-1.62, -3.49, -1.95, -3.83, -2.58, -4.44, # 5
					-1.62, -3.70, -1.95, -4.04, -2.58, -4.67, # 6
					-1.62, -3.90, -1.95, -4.23, -2.58, -4.88, # 7
					-1.62, -4.09, -1.95, -4.43, -2.58, -5.07, # 8
					-1.62, -4.26, -1.95, -4.61, -2.58, -5.25, # 9
					-1.62, -4.42, -1.95, -4.76, -2.58, -5.44), ncol = 6, byrow = T)	 # 10
				}
			}
		# Case II: restricted intercept, no trend
		else if (case == 2)	{
			#		0.10         0.05		   0.010
			# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
			fmat <- matrix(c(3.980, 3.980, 4.945, 4.945, 7.350, 7.350, # 0
				3.223, 3.757, 3.957, 4.530, 5.763, 6.480, # 1
				2.845, 3.623, 3.478, 4.335, 4.948, 6.028, # 2
				2.618, 3.532, 3.164, 4.194, 4.428, 5.816, # 3
				2.460, 3.460, 2.947, 4.088, 4.093, 5.532, # 4
				2.331, 3.417, 2.804, 4.013, 3.900, 5.419, # 5
				2.254, 3.388, 2.685, 3.960, 3.713, 5.326, # 6
				2.196, 3.370, 2.597, 3.907, 3.599, 5.230), ncol = 6, byrow = T) # 7
			if (!is.null(tstat)) {
				tnote <- "Critical values do not currently exist for Case II."
				}
			}
		# Case III: unrestricted intercept, no trend
		else if (case == 3)	{
			#		0.10         0.05		   0.010
			# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
			fmat <- matrix(c(6.810, 6.810, 8.640, 8.640, 13.290, 13.290, # 0
				4.225, 5.050, 5.290, 6.175, 7.870, 8.960, # 1
				3.393, 4.410, 4.183, 5.333, 6.140, 7.607, # 2
				2.958, 4.100, 3.615, 4.913, 5.198, 6.845, # 3
				2.696, 3.898, 3.276, 4.630, 4.590, 6.368, # 4
				2.508, 3.763, 3.037, 4.443, 4.257, 6.040, # 5
				2.387, 3.671, 2.864, 4.324, 4.016, 5.797, # 6
				2.300, 3.606, 2.753, 4.209, 3.841, 5.686), ncol = 6, byrow = T) # 7
			if (!is.null(tstat)) {
				tnote <- "Small-sample critical values not provided for Case III. Asymptotic critical values used."
				#		0.10         0.05		   0.010
				# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
				tmat <- matrix(c(-2.57, -2.57, -2.86, -2.86, -3.43, -3.43, # 0
					-2.57, -2.91, -2.86, -3.22, -3.43, -3.82, # 1
					-2.57, -3.21, -2.86, -3.53, -3.43, -4.10, # 2
					-2.57, -3.46, -2.86, -3.78, -3.43, -4.37, # 3
					-2.57, -3.66, -2.86, -3.99, -3.43, -4.60, # 4
					-2.57, -3.86, -2.86, -4.19, -3.43, -4.79, # 5
					-2.57, -4.04, -2.86, -4.38, -3.43, -4.99, # 6
					-2.57, -4.23, -2.86, -4.57, -3.43, -5.19, # 7
					-2.57, -4.40, -2.86, -4.72, -3.43, -5.37, # 8
					-2.57, -4.56, -2.86, -4.88, -3.42, -5.54, # 9
					-2.57, -4.69, -2.86, -5.03, -3.43, -5.68), ncol = 6, byrow = T) # 10
				}
			}
		# Case IV: unrestricted intercept, restricted trend
		else if (case == 4)	{
			#		0.10         0.05		   0.010
			# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
			fmat <- matrix(c(5.690, 5.690, 6.900, 6.900, 9.975, 9.975, # 0
				4.380, 4.867, 5.233, 5.777, 7.477, 8.213, # 1
				3.698, 4.420, 4.433, 5.245, 6.328, 7.408, # 2
				3.290, 4.176, 3.936, 4.918, 5.654, 6.926, # 3
				3.035, 3.997, 3.578, 4.668, 5.147, 6.617, # 4
				2.831, 3.879, 3.353, 4.500, 4.849, 6.511, # 5
				2.685, 3.785, 3.174, 4.383, 4.629, 5.698, # 6
				2.578, 3.710, 3.057, 4.319, 4.489, 5.064), ncol = 6, byrow = T)	 # 7
			if (!is.null(tstat)) {
				tnote <- "Critical values do not currently exist for Case IV."
				}
			}
		# Case V: unrestricted intercept, unrestricted trend
		else { # case == 5
			#		0.10         0.05		   0.010
			# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
			fmat <- matrix(c(10.240, 10.240, 12.580, 12.580, 18.020, 18.020, # 0
				5.950, 6.680, 7.210, 8.055, 10.365, 11.295, # 1
				4.517, 5.480, 5.457, 6.570, 7.643, 9.063, # 2
				3.800, 4.888, 4.568, 5.795, 6.380, 7.730, # 3
				3.374, 4.512, 4.036, 5.304, 5.604, 7.172, # 4
				3.087, 4.277, 3.673, 5.002, 5.095, 6.770, # 5
				2.879, 4.114, 3.426, 4.790, 4.704, 6.537, # 6
				2.729, 3.985, 3.251, 4.640, 4.459, 6.206), ncol = 6, byrow = T) # 7
			if (!is.null(tstat)) {
				tnote <- "Small-sample critical values not provided for Case V. Asymptotic critical values used."
				#		0.10         0.05		   0.010
				# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
				tmat <- matrix(c(-3.13, -3.13, -3.41, -3.41, -3.96, -3.97, # 0
					-3.13, -3.40, -3.41, -3.69, -3.96, -4.26, # 1
					-3.13, -3.63, -3.41, -3.95, -3.96, -4.53, # 2
					-3.13, -3.84, -3.41, -4.16, -3.96, -4.73, # 3
					-3.13, -4.04, -3.41, -4.36, -3.96, -4.96, # 4
					-3.13, -4.21, -3.41, -4.52, -3.96, -5.13, # 5
					-3.13, -4.37, -3.41, -4.69, -3.96, -5.31, # 6
					-3.13, -4.53, -3.41, -4.85, -3.96, -5.49, # 7
					-3.13, -4.68, -3.41, -5.01, -3.96, -5.65, # 8
					-3.13, -4.82, -3.41, -5.15, -3.96, -5.79, # 9
					-3.13, -4.96, -3.41, -5.29, -3.96, -5.94), ncol = 6, byrow = T) # 10
				}
			}
		}
	else if (obs <= 40)	{
		# Case I: no intercept, no trend
		if (case == 1) {
			# Doesn't exist! Use asym. crit values
			fnote <- "Small-sample critical values are not provided for Case I. Asymptotic critical values used."
			#		0.10        0.05	   0.010
			# 	 I(0)   I(1)  I(0)  I(1)  I(0) I(1)
			fmat <- matrix(c(3.00, 3.00, 4.20, 4.20, 7.17, 7.17, # 0
				2.44, 3.28, 3.15, 4.11, 4.81, 6.02, # 1
				2.17, 3.19, 2.72, 3.83, 3.88, 5.30, # 2
				2.01, 3.10, 2.45, 3.63, 3.42, 4.84, # 3
				1.90, 3.01, 2.26, 3.48, 3.07, 4.44, # 4
				1.81, 2.93, 2.14, 3.34, 2.82, 4.21, # 5
				1.75, 2.87, 2.04, 3.24, 2.66, 4.05, # 6
				1.70, 2.83, 1.97, 3.18, 2.54, 3.91, # 7
				1.66, 2.79, 1.91, 3.11, 2.45, 3.79, # 8
				1.63, 2.75, 1.86, 3.05, 2.34, 3.68, # 9
				1.60, 2.72, 1.82, 2.99, 2.26, 3.60), ncol = 6, byrow = T) # 10
			if (!is.null(tstat)) {
				tnote <- "Small-sample critical values not provided for Case I. Asymptotic critical values used."
				#		0.10         0.05		   0.010
				# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
				tmat <- matrix(c(-1.62, -1.62, -1.95, -1.95, -2.58, -2.58, # 0
					-1.62, -2.28, -1.95, -2.60, -2.58, -3.22, # 1
					-1.62, -2.68, -1.95, -3.02, -2.58, -3.66, # 2
					-1.62, -3.00, -1.95, -3.33, -2.58, -3.97, # 3
					-1.62, -3.26, -1.95, -3.60, -2.58, -4.23, # 4
					-1.62, -3.49, -1.95, -3.83, -2.58, -4.44, # 5
					-1.62, -3.70, -1.95, -4.04, -2.58, -4.67, # 6
					-1.62, -3.90, -1.95, -4.23, -2.58, -4.88, # 7
					-1.62, -4.09, -1.95, -4.43, -2.58, -5.07, # 8
					-1.62, -4.26, -1.95, -4.61, -2.58, -5.25, # 9
					-1.62, -4.42, -1.95, -4.76, -2.58, -5.44), ncol = 6, byrow = T) # 10
				}
			}
		# Case II: restricted intercept, no trend
		else if (case == 2)	{
			#		0.10         0.05		   0.010
			# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
			fmat <- matrix(c(3.955, 3.955, 4.960, 4.960, 7.220, 7.220, # 0
				3.210, 3.730, 3.937, 4.523, 5.593, 6.333, # 1
				2.835, 3.585, 3.435, 4.260, 4.770, 5.855, # 2
				2.592, 3.454, 3.100, 4.088, 4.310, 5.544, # 3
				2.427, 3.395, 2.893, 4.000, 3.967, 5.455, # 4
				2.306, 3.353, 2.734, 3.920, 3.657, 5.256, # 5
				2.218, 3.314, 2.618, 3.863, 3.505, 5.121, # 6
				2.152, 3.296, 2.523, 3.829, 3.402, 5.031), ncol = 6, byrow = T) # 7
			if (!is.null(tstat)) {
				tnote <- "Critical values do not currently exist for Case II."
				}
			}
		# Case III: unrestricted intercept, no trend
		else if (case == 3)	{
			#		0.10         0.05		   0.010
			# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
			fmat <- matrix(c(6.760, 6.760, 8.570, 8.570, 13.070, 13.070, # 0
				4.235, 5.000, 5.260, 6.160, 7.625, 8.825, # 1
				3.373, 4.377, 4.133, 5.260, 5.893, 7.337, # 2
				2.933, 4.020, 3.548, 4.803, 5.018, 6.610, # 3
				2.660, 3.838, 3.202, 4.544, 4.428, 6.250, # 4
				2.483, 3.708, 2.962, 4.338, 4.045, 5.898, # 5
				2.353, 3.599, 2.797, 4.211, 3.800, 5.643, # 6
				2.260, 3.534, 2.676, 4.130, 3.644, 5.464), ncol = 6, byrow = T) # 7
			if (!is.null(tstat)) {
				tnote <- "Small-sample critical values not provided for Case III. Asymptotic critical values used."
				#		0.10         0.05		   0.010
				# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
				tmat <- matrix(c(-2.57, -2.57, -2.86, -2.86, -3.43, -3.43, # 0
					-2.57, -2.91, -2.86, -3.22, -3.43, -3.82, # 1
					-2.57, -3.21, -2.86, -3.53, -3.43, -4.10, # 2
					-2.57, -3.46, -2.86, -3.78, -3.43, -4.37, # 3
					-2.57, -3.66, -2.86, -3.99, -3.43, -4.60, # 4
					-2.57, -3.86, -2.86, -4.19, -3.43, -4.79, # 5
					-2.57, -4.04, -2.86, -4.38, -3.43, -4.99, # 6
					-2.57, -4.23, -2.86, -4.57, -3.43, -5.19, # 7
					-2.57, -4.40, -2.86, -4.72, -3.43, -5.37, # 8
					-2.57, -4.56, -2.86, -4.88, -3.42, -5.54, # 9
					-2.57, -4.69, -2.86, -5.03, -3.43, -5.68), ncol = 6, byrow = T) # 10
				}
			}
		# Case IV: unrestricted intercept, restricted trend
		else if (case == 4)	{
			#		0.10         0.05		   0.010
			# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
			fmat <- matrix(c(5.680, 5.680, 6.870, 6.870, 9.575, 9.575, # 0
				4.343, 4.823, 5.180, 5.733, 7.207, 7.860, # 1
				3.663, 4.378, 4.360, 5.138, 5.980, 6.973, # 2
				3.264, 4.094, 3.850, 4.782, 5.258, 6.526, # 3
				2.985, 3.918, 3.512, 4.587, 4.763, 6.200, # 4
				2.781, 3.813, 3.257, 4.431, 4.427, 5.837, # 5
				2.634, 3.719, 3.070, 4.309, 4.154, 5.699, # 6
				2.517, 3.650, 2.933, 4.224, 3.971, 5.486), ncol = 6, byrow = T) # 7
			if (!is.null(tstat)) {
				tnote <- "Critical values do not currently exist for Case IV."
				}
			}
		# Case V: unrestricted intercept, unrestricted trend
		else { # case == 5
			#		0.10         0.05		   0.010
			# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
			fmat <- matrix(c(10.160, 10.160, 12.510, 12.510, 17.910, 17.910, # 0
				5.915, 6.630, 7.135, 7.980, 10.150, 11.230, # 1
				4.477, 5.420, 5.387, 6.437, 7.527, 8.803, # 2
				3.760, 4.795, 4.510, 5.643, 6.238, 7.740, # 3
				3.334, 4.438, 3.958, 5.226, 5.376, 7.092, # 4
				3.032, 4.213, 3.577, 4.923, 4.885, 6.550, # 5
				2.831, 4.040, 3.327, 4.700, 4.527, 6.263, # 6
				2.668, 3.920, 3.121, 4.564, 4.310, 5.965), ncol = 6, byrow = T) # 7
			if (!is.null(tstat)) {
				tnote <- "Small-sample critical values not provided for Case V. Asymptotic critical values used."
				#		0.10         0.05		   0.010
				# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
				tmat <- matrix(c(-3.13, -3.13, -3.41, -3.41, -3.96, -3.97, # 0
					-3.13, -3.40, -3.41, -3.69, -3.96, -4.26, # 1
					-3.13, -3.63, -3.41, -3.95, -3.96, -4.53, # 2
					-3.13, -3.84, -3.41, -4.16, -3.96, -4.73, # 3
					-3.13, -4.04, -3.41, -4.36, -3.96, -4.96, # 4
					-3.13, -4.21, -3.41, -4.52, -3.96, -5.13, # 5
					-3.13, -4.37, -3.41, -4.69, -3.96, -5.31, # 6
					-3.13, -4.53, -3.41, -4.85, -3.96, -5.49, # 7
					-3.13, -4.68, -3.41, -5.01, -3.96, -5.65, # 8
					-3.13, -4.82, -3.41, -5.15, -3.96, -5.79, # 9
					-3.13, -4.96, -3.41, -5.29, -3.96, -5.94), ncol = 6, byrow = T) # 10
				}
			}
		}
	else if (obs <= 45) {
		# Case I: no intercept, no trend
		if (case == 1)	{
			# Doesn't exist! Use asym. crit values
			fnote <- "Small-sample critical values are not provided for Case I. Asymptotic critical values used."
			#		0.10        0.05	   0.010
			# 	 I(0)   I(1)  I(0)  I(1)  I(0) I(1)
			fmat <- matrix(c(3.00, 3.00, 4.20, 4.20, 7.17, 7.17, # 0
				2.44, 3.28, 3.15, 4.11, 4.81, 6.02, # 1
				2.17, 3.19, 2.72, 3.83, 3.88, 5.30, # 2
				2.01, 3.10, 2.45, 3.63, 3.42, 4.84, # 3
				1.90, 3.01, 2.26, 3.48, 3.07, 4.44, # 4
				1.81, 2.93, 2.14, 3.34, 2.82, 4.21, # 5
				1.75, 2.87, 2.04, 3.24, 2.66, 4.05, # 6
				1.70, 2.83, 1.97, 3.18, 2.54, 3.91, # 7
				1.66, 2.79, 1.91, 3.11, 2.45, 3.79, # 8
				1.63, 2.75, 1.86, 3.05, 2.34, 3.68, # 9
				1.60, 2.72, 1.82, 2.99, 2.26, 3.60), ncol = 6, byrow = T) # 10
			if (!is.null(tstat)) {
				tnote <- "Small-sample critical values not provided for Case I. Asymptotic critical values used."
				#		0.10         0.05		   0.010
				# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
				tmat <- matrix(c(-1.62, -1.62, -1.95, -1.95, -2.58, -2.58, # 0
					-1.62, -2.28, -1.95, -2.60, -2.58, -3.22, # 1
					-1.62, -2.68, -1.95, -3.02, -2.58, -3.66, # 2
					-1.62, -3.00, -1.95, -3.33, -2.58, -3.97, # 3
					-1.62, -3.26, -1.95, -3.60, -2.58, -4.23, # 4
					-1.62, -3.49, -1.95, -3.83, -2.58, -4.44, # 5
					-1.62, -3.70, -1.95, -4.04, -2.58, -4.67, # 6
					-1.62, -3.90, -1.95, -4.23, -2.58, -4.88, # 7
					-1.62, -4.09, -1.95, -4.43, -2.58, -5.07, # 8
					-1.62, -4.26, -1.95, -4.61, -2.58, -5.25, # 9
					-1.62, -4.42, -1.95, -4.76, -2.58, -5.44), ncol = 6, byrow = T) # 10
				}
			}
		# Case II: restricted intercept, no trend
		else if (case == 2)	{
			#		0.10         0.05		   0.010
			# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
			fmat <- matrix(c(3.950, 3.950, 4.895, 4.895, 7.265, 7.265, # 0
				3.190, 3.730, 3.877, 4.460, 5.607, 6.193, # 1
				2.788, 3.540, 3.368, 4.203, 4.800, 5.725, # 2
				2.560, 3.428, 3.078, 4.022, 4.270, 5.412, # 3
				2.402, 3.345, 2.850, 3.905, 3.892, 5.173, # 4
				2.276, 3.297, 2.694, 3.829, 3.674, 5.019, # 5
				2.188, 3.254, 2.591, 3.766, 3.540, 4.931, # 6
				2.131, 3.223, 2.504, 3.723, 3.383, 4.832), ncol = 6, byrow = T) # 7
			if (!is.null(tstat)) {
				tnote <- "Critical values do not currently exist for Case II."
				}
			}
		# Case III: unrestricted intercept, no trend
		else if (case == 3)	{
			#		0.10         0.05		   0.010
			# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
			fmat <- matrix(c(6.760, 6.760, 8.590, 8.590, 12.930, 12.930, # 0
				4.225, 5.020, 5.235, 6.135, 7.740, 8.650, # 1
				3.330, 4.347, 4.083, 5.207, 5.920, 7.197, # 2
				2.893, 3.983, 3.535, 4.733, 4.983, 6.423, # 3
				2.638, 3.772, 3.178, 4.450, 4.394, 5.914, # 4
				2.458, 3.647, 2.922, 4.268, 4.030, 5.598, # 5
				2.327, 3.541, 2.764, 4.123, 3.790, 5.411, # 6
				2.238, 3.461, 2.643, 4.004, 3.595, 5.225), ncol = 6, byrow = T) # 7
			if (!is.null(tstat)) {
				tnote <- "Small-sample critical values not provided for Case III. Asymptotic critical values used."
				#		0.10         0.05		   0.010
				# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
				tmat <- matrix(c(-2.57, -2.57, -2.86, -2.86, -3.43, -3.43, # 0
					-2.57, -2.91, -2.86, -3.22, -3.43, -3.82, # 1
					-2.57, -3.21, -2.86, -3.53, -3.43, -4.10, # 2
					-2.57, -3.46, -2.86, -3.78, -3.43, -4.37, # 3
					-2.57, -3.66, -2.86, -3.99, -3.43, -4.60, # 4
					-2.57, -3.86, -2.86, -4.19, -3.43, -4.79, # 5
					-2.57, -4.04, -2.86, -4.38, -3.43, -4.99, # 6
					-2.57, -4.23, -2.86, -4.57, -3.43, -5.19, # 7
					-2.57, -4.40, -2.86, -4.72, -3.43, -5.37, # 8
					-2.57, -4.56, -2.86, -4.88, -3.42, -5.54, # 9
					-2.57, -4.69, -2.86, -5.03, -3.43, -5.68), ncol = 6, byrow = T) # 10
				}
			}
		# Case IV: unrestricted intercept, restricted trend
		else if (case == 4)	{
			#		0.10         0.05		   0.010
			# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
			fmat <- matrix(c(5.625, 5.625, 6.750, 6.750, 9.555, 9.555, # 0
				4.300, 4.780, 5.130, 5.680, 7.133, 7.820, # 1
				3.625, 4.330, 4.335, 5.078, 5.878, 6.870, # 2
				3.226, 4.054, 3.822, 4.714, 5.150, 6.280, # 3
				2.950, 3.862, 3.470, 4.470, 4.628, 5.865, # 4
				2.750, 3.739, 3.211, 4.309, 4.251, 5.596, # 5
				2.606, 3.644, 3.025, 4.198, 3.998, 5.463, # 6
				2.484, 3.570, 2.899, 4.087, 3.829, 5.313), ncol = 6, byrow = T) # 7
			if (!is.null(tstat)) {
				tnote <- "Critical values do not currently exist for Case IV."
				}
			}
		# Case V: unrestricted intercept, unrestricted trend
		else { # case == 5
			#		0.10         0.05		   0.010
			# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
			fmat <- matrix(c(10.150, 10.150, 12.400, 12.400, 17.500, 17.500, # 0
				5.880, 6.640, 7.080, 7.910, 9.890, 10.965, # 1
				4.437, 5.377, 5.360, 6.373, 7.317, 8.720, # 2
				3.740, 4.780, 4.450, 5.560, 6.053, 7.458, # 3
				3.298, 4.378, 3.890, 5.104, 5.224, 6.696, # 4
				3.012, 4.147, 3.532, 4.800, 4.715, 6.293, # 5
				2.796, 3.970, 3.267, 4.584, 4.364, 6.006, # 6
				2.635, 3.838, 3.091, 4.413, 4.109, 5.785), ncol = 6, byrow = T) # 7
			if (!is.null(tstat)) {
				tnote <- "Small-sample critical values not provided for Case V. Asymptotic critical values used."
				#		0.10         0.05		   0.010
				# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
				tmat <- matrix(c(-3.13, -3.13, -3.41, -3.41, -3.96, -3.97, # 0
					-3.13, -3.40, -3.41, -3.69, -3.96, -4.26, # 1
					-3.13, -3.63, -3.41, -3.95, -3.96, -4.53, # 2
					-3.13, -3.84, -3.41, -4.16, -3.96, -4.73, # 3
					-3.13, -4.04, -3.41, -4.36, -3.96, -4.96, # 4
					-3.13, -4.21, -3.41, -4.52, -3.96, -5.13, # 5
					-3.13, -4.37, -3.41, -4.69, -3.96, -5.31, # 6
					-3.13, -4.53, -3.41, -4.85, -3.96, -5.49, # 7
					-3.13, -4.68, -3.41, -5.01, -3.96, -5.65, # 8
					-3.13, -4.82, -3.41, -5.15, -3.96, -5.79, # 9
					-3.13, -4.96, -3.41, -5.29, -3.96, -5.94), ncol = 6, byrow = T) # 10
				}
			}
		}
	else if (obs <= 50)	{
		# Case I: no intercept, no trend
		if (case == 1) {
			# Doesn't exist! Use asym. crit values
			fnote <- "Small-sample critical values are not provided for Case I. Asymptotic critical values used."
			#		0.10        0.05	   0.010
			# 	 I(0)   I(1)  I(0)  I(1)  I(0) I(1)
			fmat <- matrix(c(3.00, 3.00, 4.20, 4.20, 7.17, 7.17, # 0
				2.44, 3.28, 3.15, 4.11, 4.81, 6.02, # 1
				2.17, 3.19, 2.72, 3.83, 3.88, 5.30, # 2
				2.01, 3.10, 2.45, 3.63, 3.42, 4.84, # 3
				1.90, 3.01, 2.26, 3.48, 3.07, 4.44, # 4
				1.81, 2.93, 2.14, 3.34, 2.82, 4.21, # 5
				1.75, 2.87, 2.04, 3.24, 2.66, 4.05, # 6
				1.70, 2.83, 1.97, 3.18, 2.54, 3.91, # 7
				1.66, 2.79, 1.91, 3.11, 2.45, 3.79, # 8
				1.63, 2.75, 1.86, 3.05, 2.34, 3.68, # 9
				1.60, 2.72, 1.82, 2.99, 2.26, 3.60), ncol = 6, byrow = T) # 10
			if (!is.null(tstat)) {
				tnote <- "Small-sample critical values not provided for Case I. Asymptotic critical values used."
				#		0.10         0.05		   0.010
				# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
				tmat <- matrix(c(-1.62, -1.62, -1.95, -1.95, -2.58, -2.58, # 0
					-1.62, -2.28, -1.95, -2.60, -2.58, -3.22, # 1
					-1.62, -2.68, -1.95, -3.02, -2.58, -3.66, # 2
					-1.62, -3.00, -1.95, -3.33, -2.58, -3.97, # 3
					-1.62, -3.26, -1.95, -3.60, -2.58, -4.23, # 4
					-1.62, -3.49, -1.95, -3.83, -2.58, -4.44, # 5
					-1.62, -3.70, -1.95, -4.04, -2.58, -4.67, # 6
					-1.62, -3.90, -1.95, -4.23, -2.58, -4.88, # 7
					-1.62, -4.09, -1.95, -4.43, -2.58, -5.07, # 8
					-1.62, -4.26, -1.95, -4.61, -2.58, -5.25, # 9
					-1.62, -4.42, -1.95, -4.76, -2.58, -5.44), ncol = 6, byrow = T) # 10
				}
			}
		# Case II: restricted intercept, no trend
		else if (case == 2)	{
			#		0.10         0.05		   0.010
			# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
			fmat <- matrix(c(3.935, 3.935, 4.815, 4.815, 7.065, 7.065, # 0
				3.177, 3.653, 3.860, 4.440, 5.503, 6.240, # 1
				2.788, 3.513, 3.368, 4.178, 4.695, 5.758, # 2
				2.538, 3.398, 3.048, 4.002, 4.188, 5.328, # 3
				2.372, 3.320, 2.823, 3.872, 3.845, 5.150, # 4
				2.259, 3.264, 2.670, 3.781, 3.593, 4.981, # 5
				2.170, 3.220, 2.550, 3.708, 3.424, 4.880, # 6
				2.099, 3.181, 2.457, 3.650, 3.282, 4.730), ncol = 6, byrow = T) # 7
			if (!is.null(tstat)) {
				tnote <- "Critical values do not currently exist for Case II."
				}
			}
		# Case III: unrestricted intercept, no trend
		else if (case == 3)	{
			#		0.10         0.05		   0.010
			# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
			fmat <- matrix(c(6.740, 6.740, 8.510, 8.510, 12.730, 12.730, # 0
				4.190, 4.940, 5.220, 6.070, 7.560, 8.685, # 1
				3.333, 4.313, 4.070, 5.190, 5.817, 7.303, # 2
				2.873, 3.973, 3.500, 4.700, 4.865, 6.360, # 3
				2.614, 3.746, 3.136, 4.416, 4.306, 5.874, # 4
				2.435, 3.600, 2.900, 4.218, 3.955, 5.583, # 5
				2.309, 3.507, 2.726, 4.057, 3.656, 5.331, # 6
				2.205, 3.421, 2.593, 3.941, 3.498, 5.149), ncol = 6, byrow = T) # 7
			if (!is.null(tstat)) {
				tnote <- "Small-sample critical values not provided for Case III. Asymptotic critical values used."
				#		0.10         0.05		   0.010
				# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
				tmat <- matrix(c(-2.57, -2.57, -2.86, -2.86, -3.43, -3.43, # 0
					-2.57, -2.91, -2.86, -3.22, -3.43, -3.82, # 1
					-2.57, -3.21, -2.86, -3.53, -3.43, -4.10, # 2
					-2.57, -3.46, -2.86, -3.78, -3.43, -4.37, # 3
					-2.57, -3.66, -2.86, -3.99, -3.43, -4.60, # 4
					-2.57, -3.86, -2.86, -4.19, -3.43, -4.79, # 5
					-2.57, -4.04, -2.86, -4.38, -3.43, -4.99, # 6
					-2.57, -4.23, -2.86, -4.57, -3.43, -5.19, # 7
					-2.57, -4.40, -2.86, -4.72, -3.43, -5.37, # 8
					-2.57, -4.56, -2.86, -4.88, -3.42, -5.54, # 9
					-2.57, -4.69, -2.86, -5.03, -3.43, -5.68), ncol = 6, byrow = T) # 10
				}
			}
		# Case IV: unrestricted intercept, restricted trend
		else if (case == 4)	{
			#		0.10         0.05		   0.010
			# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
			fmat <- matrix(c(5.570, 5.570, 6.685, 6.685, 9.320, 9.320, # 0
				4.230, 4.740, 5.043, 5.607, 7.017, 7.727, # 1
				3.573, 4.288, 4.225, 5.030, 5.805, 6.790, # 2
				3.174, 4.004, 3.730, 4.666, 5.050, 6.182, # 3
				2.905, 3.822, 3.383, 4.432, 4.557, 5.793, # 4
				2.703, 3.697, 3.149, 4.293, 4.214, 5.520, # 5
				2.550, 3.609, 2.975, 4.143, 3.983, 5.345, # 6
				2.440, 3.523, 2.832, 4.012, 3.762, 5.172), ncol = 6, byrow = T) # 7
			if (!is.null(tstat)) {
				tnote <- "Critical values do not currently exist for Case IV."
				}
			}
		# Case V: unrestricted intercept, unrestricted trend
		else { # case == 5
			#		0.10         0.05		   0.010
			# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
			fmat <- matrix(c(10.020, 10.020, 12.170, 12.170, 17.530, 17.530, # 0
				5.780, 6.540, 6.985, 7.860, 9.895, 10.965, # 1
				4.380, 5.350, 5.247, 6.303, 7.337, 8.643, # 2
				3.673, 4.715, 4.368, 5.545, 5.995, 7.335, # 3
				3.240, 4.350, 3.834, 5.064, 5.184, 6.684, # 4
				2.950, 4.110, 3.480, 4.782, 4.672, 6.232, # 5
				2.750, 3.944, 3.229, 4.536, 4.310, 5.881, # 6
				2.590, 3.789, 3.039, 4.339, 4.055, 5.640), ncol = 6, byrow = T) # 7
			if (!is.null(tstat)) {
				tnote <- "Small-sample critical values not provided for Case V. Asymptotic critical values used."
				#		0.10         0.05		   0.010
				# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
				tmat <- matrix(c(-3.13, -3.13, -3.41, -3.41, -3.96, -3.97, # 0
					-3.13, -3.40, -3.41, -3.69, -3.96, -4.26, # 1
					-3.13, -3.63, -3.41, -3.95, -3.96, -4.53, # 2
					-3.13, -3.84, -3.41, -4.16, -3.96, -4.73, # 3
					-3.13, -4.04, -3.41, -4.36, -3.96, -4.96, # 4
					-3.13, -4.21, -3.41, -4.52, -3.96, -5.13, # 5
					-3.13, -4.37, -3.41, -4.69, -3.96, -5.31, # 6
					-3.13, -4.53, -3.41, -4.85, -3.96, -5.49, # 7
					-3.13, -4.68, -3.41, -5.01, -3.96, -5.65, # 8
					-3.13, -4.82, -3.41, -5.15, -3.96, -5.79, # 9
					-3.13, -4.96, -3.41, -5.29, -3.96, -5.94), ncol = 6, byrow = T) # 10
				}
			}
		}
	else if (obs <= 55)	{
		# Case I: no intercept, no trend
		if (case == 1)	{
			# Doesn't exist! Use asym. crit values
			fnote <- "Small-sample critical values are not provided for Case I. Asymptotic critical values used."
			#		0.10        0.05	   0.010
			# 	 I(0)   I(1)  I(0)  I(1)  I(0) I(1)
			fmat <- matrix(c(3.00, 3.00, 4.20, 4.20, 7.17, 7.17, # 0
				2.44, 3.28, 3.15, 4.11, 4.81, 6.02, # 1
				2.17, 3.19, 2.72, 3.83, 3.88, 5.30, # 2
				2.01, 3.10, 2.45, 3.63, 3.42, 4.84, # 3
				1.90, 3.01, 2.26, 3.48, 3.07, 4.44, # 4
				1.81, 2.93, 2.14, 3.34, 2.82, 4.21, # 5
				1.75, 2.87, 2.04, 3.24, 2.66, 4.05, # 6
				1.70, 2.83, 1.97, 3.18, 2.54, 3.91, # 7
				1.66, 2.79, 1.91, 3.11, 2.45, 3.79, # 8
				1.63, 2.75, 1.86, 3.05, 2.34, 3.68, # 9
				1.60, 2.72, 1.82, 2.99, 2.26, 3.60), ncol = 6, byrow = T) # 10
			if (!is.null(tstat)) {
				tnote <- "Small-sample critical values not provided for Case I. Asymptotic critical values used."
				#		0.10         0.05		   0.010
				# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
				tmat <- matrix(c(-1.62, -1.62, -1.95, -1.95, -2.58, -2.58, # 0
					-1.62, -2.28, -1.95, -2.60, -2.58, -3.22, # 1
					-1.62, -2.68, -1.95, -3.02, -2.58, -3.66, # 2
					-1.62, -3.00, -1.95, -3.33, -2.58, -3.97, # 3
					-1.62, -3.26, -1.95, -3.60, -2.58, -4.23, # 4
					-1.62, -3.49, -1.95, -3.83, -2.58, -4.44, # 5
					-1.62, -3.70, -1.95, -4.04, -2.58, -4.67, # 6
					-1.62, -3.90, -1.95, -4.23, -2.58, -4.88, # 7
					-1.62, -4.09, -1.95, -4.43, -2.58, -5.07, # 8
					-1.62, -4.26, -1.95, -4.61, -2.58, -5.25, # 9
					-1.62, -4.42, -1.95, -4.76, -2.58, -5.44), ncol = 6, byrow = T) # 10
				}
			}
		# Case II: restricted intercept, no trend
		else if (case == 2)	{
			#		0.10         0.05		   0.010
			# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
			fmat <- matrix(c(3.900, 3.900, 4.795, 4.795, 6.965, 6.965, # 0
				3.143, 3.670, 3.790, 4.393, 5.377, 6.047, # 1
				2.748, 3.495, 3.303, 4.100, 4.610, 5.563, # 2
				2.508, 3.356, 2.982, 3.942, 4.118, 5.200, # 3
				2.345, 3.280, 2.763, 3.813, 3.738, 4.947, # 4
				2.226, 3.241, 2.617, 3.743, 3.543, 4.839, # 5
				2.139, 3.204, 2.490, 3.658, 3.330, 4.708, # 6
				2.069, 3.148, 2.414, 3.608, 3.194, 4.562), ncol = 6, byrow = T) # 7
			if (!is.null(tstat)) {
				tnote <- "Critical values do not currently exist for Case II."
				}
			}
		# Case III: unrestricted intercept, no trend
		else if (case == 3)	{
			#		0.10         0.05		   0.010
			# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
			fmat <- matrix(c(6.700, 6.700, 8.390, 8.390, 12.700, 12.700, # 0
				4.155, 4.925, 5.125, 6.045, 7.435, 8.460, # 1
				3.280, 4.273, 3.987, 5.090, 5.707, 6.977, # 2
				2.843, 3.920, 3.408, 4.623, 4.828, 6.195, # 3
				2.578, 3.710, 3.068, 4.334, 4.244, 5.726, # 4
				2.393, 3.583, 2.848, 4.160, 3.928, 5.408, # 5
				2.270, 3.486, 2.676, 3.999, 3.636, 5.169, # 6
				2.181, 3.398, 2.556, 3.904, 3.424, 4.989), ncol = 6, byrow = T) # 7
			if (!is.null(tstat)) {
				tnote <- "Small-sample critical values not provided for Case III. Asymptotic critical values used."
				#		0.10         0.05		   0.010
				# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
				tmat <- matrix(c(-2.57, -2.57, -2.86, -2.86, -3.43, -3.43, # 0
					-2.57, -2.91, -2.86, -3.22, -3.43, -3.82, # 1
					-2.57, -3.21, -2.86, -3.53, -3.43, -4.10, # 2
					-2.57, -3.46, -2.86, -3.78, -3.43, -4.37, # 3
					-2.57, -3.66, -2.86, -3.99, -3.43, -4.60, # 4
					-2.57, -3.86, -2.86, -4.19, -3.43, -4.79, # 5
					-2.57, -4.04, -2.86, -4.38, -3.43, -4.99, # 6
					-2.57, -4.23, -2.86, -4.57, -3.43, -5.19, # 7
					-2.57, -4.40, -2.86, -4.72, -3.43, -5.37, # 8
					-2.57, -4.56, -2.86, -4.88, -3.42, -5.54, # 9
					-2.57, -4.69, -2.86, -5.03, -3.43, -5.68), ncol = 6, byrow = T) # 10
				}
			}
		# Case IV: unrestricted intercept, restricted trend
		else if (case == 4)	{
			#		0.10         0.05		   0.010
			# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
			fmat <- matrix(c(5.570, 5.570, 6.660, 6.660, 9.300, 9.300, # 0
				4.230, 4.730, 5.013, 5.547, 6.893, 7.537, # 1
				3.553, 4.238, 4.183, 4.955, 5.678, 6.578, # 2
				3.132, 3.956, 3.692, 4.582, 4.990, 6.018, # 3
				2.868, 3.782, 3.358, 4.365, 4.455, 5.615, # 4
				2.674, 3.659, 3.131, 4.206, 4.111, 5.329, # 5
				2.538, 3.560, 2.946, 4.065, 3.870, 5.171, # 6
				2.420, 3.481, 2.791, 3.950, 3.643, 5.021), ncol = 6, byrow = T) # 7
			if (!is.null(tstat)) {
				tnote <- "Critical values do not currently exist for Case IV."
				}
			}
		# Case V: unrestricted intercept, unrestricted trend
		else { # case == 5
			#		0.10         0.05		   0.010
			# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
			fmat <- matrix(c(10.110, 10.110, 12.170, 12.170, 17.480, 17.480, # 0
				5.800, 6.515, 6.930, 7.785, 9.800, 10.675, # 1
				4.370, 5.303, 5.190, 6.223, 7.227, 8.340, # 2
				3.640, 4.670, 4.313, 5.425, 5.955, 7.225, # 3
				3.210, 4.294, 3.794, 4.986, 5.108, 6.494, # 4
				2.927, 4.068, 3.442, 4.690, 4.608, 5.977, # 5
				2.724, 3.893, 3.197, 4.460, 4.230, 5.713, # 6
				2.573, 3.760, 2.989, 4.271, 3.955, 5.474), ncol = 6, byrow = T) # 7
			if (!is.null(tstat)) {
				tnote <- "Small-sample critical values not provided for Case V. Asymptotic critical values used."
				#		0.10         0.05		   0.010
				# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
				tmat <- matrix(c(-3.13, -3.13, -3.41, -3.41, -3.96, -3.97, # 0
					-3.13, -3.40, -3.41, -3.69, -3.96, -4.26, # 1
					-3.13, -3.63, -3.41, -3.95, -3.96, -4.53, # 2
					-3.13, -3.84, -3.41, -4.16, -3.96, -4.73, # 3
					-3.13, -4.04, -3.41, -4.36, -3.96, -4.96, # 4
					-3.13, -4.21, -3.41, -4.52, -3.96, -5.13, # 5
					-3.13, -4.37, -3.41, -4.69, -3.96, -5.31, # 6
					-3.13, -4.53, -3.41, -4.85, -3.96, -5.49, # 7
					-3.13, -4.68, -3.41, -5.01, -3.96, -5.65, # 8
					-3.13, -4.82, -3.41, -5.15, -3.96, -5.79, # 9
					-3.13, -4.96, -3.41, -5.29, -3.96, -5.94), ncol = 6, byrow = T) # 10
				}
			}
		}
	else if (obs <= 60)	{
		# Case I: no intercept, no trend
		if (case == 1) {
			# Doesn't exist! Use asym. crit values
			fnote <- "Small-sample critical values are not provided for Case I. Asymptotic critical values used."
			#		0.10        0.05	   0.010
			# 	 I(0)   I(1)  I(0)  I(1)  I(0) I(1)
			fmat <- matrix(c(3.00, 3.00, 4.20, 4.20, 7.17, 7.17, # 0
				2.44, 3.28, 3.15, 4.11, 4.81, 6.02, # 1
				2.17, 3.19, 2.72, 3.83, 3.88, 5.30, # 2
				2.01, 3.10, 2.45, 3.63, 3.42, 4.84, # 3
				1.90, 3.01, 2.26, 3.48, 3.07, 4.44, # 4
				1.81, 2.93, 2.14, 3.34, 2.82, 4.21, # 5
				1.75, 2.87, 2.04, 3.24, 2.66, 4.05, # 6
				1.70, 2.83, 1.97, 3.18, 2.54, 3.91, # 7
				1.66, 2.79, 1.91, 3.11, 2.45, 3.79, # 8
				1.63, 2.75, 1.86, 3.05, 2.34, 3.68, # 9
				1.60, 2.72, 1.82, 2.99, 2.26, 3.60), ncol = 6, byrow = T) # 10
			if (!is.null(tstat)) {
				tnote <- "Small-sample critical values not provided for Case I. Asymptotic critical values used."
				#		0.10         0.05		   0.010
				# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
				tmat <- matrix(c(-1.62, -1.62, -1.95, -1.95, -2.58, -2.58, # 0
					-1.62, -2.28, -1.95, -2.60, -2.58, -3.22, # 1
					-1.62, -2.68, -1.95, -3.02, -2.58, -3.66, # 2
					-1.62, -3.00, -1.95, -3.33, -2.58, -3.97, # 3
					-1.62, -3.26, -1.95, -3.60, -2.58, -4.23, # 4
					-1.62, -3.49, -1.95, -3.83, -2.58, -4.44, # 5
					-1.62, -3.70, -1.95, -4.04, -2.58, -4.67, # 6
					-1.62, -3.90, -1.95, -4.23, -2.58, -4.88, # 7
					-1.62, -4.09, -1.95, -4.43, -2.58, -5.07, # 8
					-1.62, -4.26, -1.95, -4.61, -2.58, -5.25, # 9
					-1.62, -4.42, -1.95, -4.76, -2.58, -5.44), ncol = 6, byrow = T) # 10
				}
			}
		# Case II: restricted intercept, no trend
		else if (case == 2)	{
			#		0.10         0.05		   0.010
			# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
			fmat <- matrix(c(3.880, 3.880, 4.780, 4.780, 6.960, 6.960, # 0
				3.127, 3.650, 3.803, 4.363, 5.383, 6.033, # 1
				2.738, 3.465, 3.288, 4.070, 4.558, 5.590, # 2
				2.496, 3.346, 2.962, 3.910, 4.068, 5.250, # 3
				2.323, 3.273, 2.743, 3.792, 3.710, 4.965, # 4
				2.204, 3.210, 2.589, 3.683, 3.451, 4.764, # 5
				2.114, 3.153, 2.456, 3.598, 3.293, 4.615, # 6
				2.044, 3.104, 2.373, 3.540, 3.129, 4.507), ncol = 6, byrow = T) # 7
			if (!is.null(tstat)) {
				tnote <- "Critical values do not currently exist for Case II."
				}
			}
		# Case III: unrestricted intercept, no trend
		else if (case == 3)	{
			#		0.10         0.05		   0.010
			# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
			fmat <- matrix(c(6.700, 6.700, 8.460, 8.460, 12.490, 12.490, # 0
				4.145, 4.950, 5.125, 6.000, 7.400, 8.510, # 1
				3.270, 4.260, 4.000, 5.057, 5.697, 6.987, # 2
				2.838, 3.923, 3.415, 4.615, 4.748, 6.188, # 3
				2.568, 3.712, 3.062, 4.314, 4.176, 5.676, # 4
				2.385, 3.565, 2.817, 4.097, 3.783, 5.338, # 5
				2.253, 3.436, 2.643, 3.939, 3.531, 5.081, # 6
				2.155, 3.353, 2.513, 3.823, 3.346, 4.895), ncol = 6, byrow = T) # 7
			if (!is.null(tstat)) {
				tnote <- "Small-sample critical values not provided for Case III. Asymptotic critical values used."
				#		0.10         0.05		   0.010
				# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
				tmat <- matrix(c(-2.57, -2.57, -2.86, -2.86, -3.43, -3.43, # 0
					-2.57, -2.91, -2.86, -3.22, -3.43, -3.82, # 1
					-2.57, -3.21, -2.86, -3.53, -3.43, -4.10, # 2
					-2.57, -3.46, -2.86, -3.78, -3.43, -4.37, # 3
					-2.57, -3.66, -2.86, -3.99, -3.43, -4.60, # 4
					-2.57, -3.86, -2.86, -4.19, -3.43, -4.79, # 5
					-2.57, -4.04, -2.86, -4.38, -3.43, -4.99, # 6
					-2.57, -4.23, -2.86, -4.57, -3.43, -5.19, # 7
					-2.57, -4.40, -2.86, -4.72, -3.43, -5.37, # 8
					-2.57, -4.56, -2.86, -4.88, -3.42, -5.54, # 9
					-2.57, -4.69, -2.86, -5.03, -3.43, -5.68), ncol = 6, byrow = T) # 10
				}
			}
		# Case IV: unrestricted intercept, restricted trend
		else if (case == 4)	{
			#		0.10         0.05		   0.010
			# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
			fmat <- matrix(c(5.555, 5.555, 6.630, 6.630, 9.245, 9.245, # 0
				4.203, 4.693, 4.980, 5.527, 6.780, 7.377, # 1
				3.540, 4.235, 4.180, 4.938, 5.620, 6.503, # 2
				3.130, 3.968, 3.684, 4.584, 4.928, 5.950, # 3
				2.852, 3.773, 3.323, 4.333, 4.412, 5.545, # 4
				2.653, 3.637, 3.086, 4.154, 4.013, 5.269, # 5
				2.510, 3.519, 2.900, 3.999, 3.775, 5.086, # 6
				2.392, 3.444, 2.756, 3.892, 3.584, 4.922), ncol = 6, byrow = T) # 7
			if (!is.null(tstat)) {
				tnote <- "Critical values do not currently exist for Case IV."
				}
			}
		# Case V: unrestricted intercept, unrestricted trend
		else { # case == 5
			#		0.10         0.05		   0.010
			# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
			fmat <- matrix(c(10.030, 10.030, 12.200, 12.200, 17.020, 17.020, # 0
				5.765, 6.500, 6.905, 7.735, 9.585, 10.420, # 1
				4.350, 5.283, 5.190, 6.200, 7.057, 8.243, # 2
				3.645, 4.678, 4.298, 5.445, 5.835, 7.108, # 3
				3.200, 4.310, 3.772, 4.956, 5.066, 6.394, # 4
				2.912, 4.047, 3.407, 4.632, 4.505, 5.920, # 5
				2.709, 3.856, 3.137, 4.393, 4.117, 5.597, # 6
				2.551, 3.716, 2.956, 4.230, 3.870, 5.338), ncol = 6, byrow = T) # 7
			if (!is.null(tstat)) {
				tnote <- "Small-sample critical values not provided for Case V. Asymptotic critical values used."
				#		0.10         0.05		   0.010
				# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
				tmat <- matrix(c(-3.13, -3.13, -3.41, -3.41, -3.96, -3.97, # 0
					-3.13, -3.40, -3.41, -3.69, -3.96, -4.26, # 1
					-3.13, -3.63, -3.41, -3.95, -3.96, -4.53, # 2
					-3.13, -3.84, -3.41, -4.16, -3.96, -4.73, # 3
					-3.13, -4.04, -3.41, -4.36, -3.96, -4.96, # 4
					-3.13, -4.21, -3.41, -4.52, -3.96, -5.13, # 5
					-3.13, -4.37, -3.41, -4.69, -3.96, -5.31, # 6
					-3.13, -4.53, -3.41, -4.85, -3.96, -5.49, # 7
					-3.13, -4.68, -3.41, -5.01, -3.96, -5.65, # 8
					-3.13, -4.82, -3.41, -5.15, -3.96, -5.79, # 9
					-3.13, -4.96, -3.41, -5.29, -3.96, -5.94), ncol = 6, byrow = T) # 10
				}
			}
		}
	else if (obs <= 65)	{
		# Case I: no intercept, no trend
		if (case == 1) {
			# Doesn't exist! Use asym. crit values
			fnote <- "Small-sample critical values are not provided for Case I. Asymptotic critical values used."
			#		0.10        0.05	   0.010
			# 	 I(0)   I(1)  I(0)  I(1)  I(0) I(1)
			fmat <- matrix(c(3.00, 3.00, 4.20, 4.20, 7.17, 7.17, # 0
				2.44, 3.28, 3.15, 4.11, 4.81, 6.02, # 1
				2.17, 3.19, 2.72, 3.83, 3.88, 5.30, # 2
				2.01, 3.10, 2.45, 3.63, 3.42, 4.84, # 3
				1.90, 3.01, 2.26, 3.48, 3.07, 4.44, # 4
				1.81, 2.93, 2.14, 3.34, 2.82, 4.21, # 5
				1.75, 2.87, 2.04, 3.24, 2.66, 4.05, # 6
				1.70, 2.83, 1.97, 3.18, 2.54, 3.91, # 7
				1.66, 2.79, 1.91, 3.11, 2.45, 3.79, # 8
				1.63, 2.75, 1.86, 3.05, 2.34, 3.68, # 9
				1.60, 2.72, 1.82, 2.99, 2.26, 3.60), ncol = 6, byrow = T) # 10
			if (!is.null(tstat)) {
				tnote <- "Small-sample critical values not provided for Case I. Asymptotic critical values used."
				#		0.10         0.05		   0.010
				# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
				tmat <- matrix(c(-1.62, -1.62, -1.95, -1.95, -2.58, -2.58, # 0
					-1.62, -2.28, -1.95, -2.60, -2.58, -3.22, # 1
					-1.62, -2.68, -1.95, -3.02, -2.58, -3.66, # 2
					-1.62, -3.00, -1.95, -3.33, -2.58, -3.97, # 3
					-1.62, -3.26, -1.95, -3.60, -2.58, -4.23, # 4
					-1.62, -3.49, -1.95, -3.83, -2.58, -4.44, # 5
					-1.62, -3.70, -1.95, -4.04, -2.58, -4.67, # 6
					-1.62, -3.90, -1.95, -4.23, -2.58, -4.88, # 7
					-1.62, -4.09, -1.95, -4.43, -2.58, -5.07, # 8
					-1.62, -4.26, -1.95, -4.61, -2.58, -5.25, # 9
					-1.62, -4.42, -1.95, -4.76, -2.58, -5.44), ncol = 6, byrow = T) # 10
				}
			}
		# Case II: restricted intercept, no trend
		else if (case == 2)	{
			#		0.10         0.05		   0.010
			# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
			fmat <- matrix(c(3.880, 3.880, 4.780, 4.780, 6.825, 6.825, # 0
				3.143, 3.623, 3.787, 4.343, 5.350, 6.017, # 1
				2.740, 3.455, 3.285, 4.070, 4.538, 5.475, # 2
				2.492, 3.350, 2.976, 3.896, 4.056, 5.158, # 3
				2.335, 3.252, 2.750, 3.755, 3.725, 4.940, # 4
				2.209, 3.201, 2.596, 3.677, 3.430, 4.721, # 5
				2.120, 3.145, 2.473, 3.583, 3.225, 4.571, # 6
				2.043, 3.094, 2.373, 3.519, 3.092, 4.478), ncol = 6, byrow = T) # 7
			if (!is.null(tstat)) {
				tnote <- "Critical values do not currently exist for Case II."
				}
			}
		# Case III: unrestricted intercept, no trend
		else if (case == 3)	{
			#		0.10         0.05		   0.010
			# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
			fmat <- matrix(c(6.740, 6.740, 8.490, 8.490, 12.400, 12.400, # 0
				4.175, 4.930, 5.130, 5.980, 7.320, 8.435, # 1
				3.300, 4.250, 4.010, 5.080, 5.583, 6.853, # 2
				2.843, 3.923, 3.435, 4.583, 4.690, 6.143, # 3
				2.574, 3.682, 3.068, 4.274, 4.188, 5.694, # 4
				2.397, 3.543, 2.835, 4.090, 3.783, 5.300, # 5
				2.256, 3.430, 2.647, 3.921, 3.501, 5.051, # 6
				2.156, 3.334, 2.525, 3.808, 3.310, 4.871), ncol = 6, byrow = T) # 7
			if (!is.null(tstat)) {
				tnote <- "Small-sample critical values not provided for Case III. Asymptotic critical values used."
				#		0.10         0.05		   0.010
				# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
				tmat <- matrix(c(-2.57, -2.57, -2.86, -2.86, -3.43, -3.43, # 0
					-2.57, -2.91, -2.86, -3.22, -3.43, -3.82, # 1
					-2.57, -3.21, -2.86, -3.53, -3.43, -4.10, # 2
					-2.57, -3.46, -2.86, -3.78, -3.43, -4.37, # 3
					-2.57, -3.66, -2.86, -3.99, -3.43, -4.60, # 4
					-2.57, -3.86, -2.86, -4.19, -3.43, -4.79, # 5
					-2.57, -4.04, -2.86, -4.38, -3.43, -4.99, # 6
					-2.57, -4.23, -2.86, -4.57, -3.43, -5.19, # 7
					-2.57, -4.40, -2.86, -4.72, -3.43, -5.37, # 8
					-2.57, -4.56, -2.86, -4.88, -3.42, -5.54, # 9
					-2.57, -4.69, -2.86, -5.03, -3.43, -5.68), ncol = 6, byrow = T) # 10
				}
			}
		# Case IV: unrestricted intercept, restricted trend
		else if (case == 4)	{
			#		0.10         0.05		   0.010
			# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
			fmat <- matrix(c(5.510, 5.510, 6.550, 6.550, 8.960, 8.960, # 0
				4.187, 4.660, 4.950, 5.467, 6.707, 7.360, # 1
				3.535, 4.208, 4.123, 4.903, 5.545, 6.453, # 2
				3.122, 3.942, 3.626, 4.538, 4.848, 5.842, # 3
				2.848, 3.743, 3.300, 4.280, 4.347, 5.552, # 4
				2.647, 3.603, 3.063, 4.123, 4.020, 5.263, # 5
				2.499, 3.490, 2.880, 3.978, 3.758, 5.040, # 6
				2.379, 3.406, 2.730, 3.879, 3.557, 4.902), ncol = 6, byrow = T) # 7
			if (!is.null(tstat)) {
				tnote <- "Critical values do not currently exist for Case IV."
				}
			}
		# Case V: unrestricted intercept, unrestricted trend
		else { # case == 5
			#		0.10         0.05		   0.010
			# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
			fmat <- matrix(c(9.970, 9.970, 11.960, 11.960, 16.850, 16.850, # 0
				5.755, 6.470, 6.890, 7.660, 9.475, 10.515, # 1
				4.353, 5.257, 5.137, 6.173, 7.013, 8.230, # 2
				3.638, 4.643, 4.268, 5.415, 5.795, 7.053, # 3
				3.196, 4.262, 3.732, 4.920, 4.974, 6.378, # 4
				2.897, 4.022, 3.372, 4.613, 4.482, 5.923, # 5
				2.690, 3.830, 3.137, 4.363, 4.111, 5.586, # 6
				2.531, 3.685, 2.924, 4.206, 3.835, 5.339), ncol = 6, byrow = T) # 7
			if (!is.null(tstat)) {
				tnote <- "Small-sample critical values not provided for Case V. Asymptotic critical values used."
				#		0.10         0.05		   0.010
				# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
				tmat <- matrix(c(-3.13, -3.13, -3.41, -3.41, -3.96, -3.97, # 0
					-3.13, -3.40, -3.41, -3.69, -3.96, -4.26, # 1
					-3.13, -3.63, -3.41, -3.95, -3.96, -4.53, # 2
					-3.13, -3.84, -3.41, -4.16, -3.96, -4.73, # 3
					-3.13, -4.04, -3.41, -4.36, -3.96, -4.96, # 4
					-3.13, -4.21, -3.41, -4.52, -3.96, -5.13, # 5
					-3.13, -4.37, -3.41, -4.69, -3.96, -5.31, # 6
					-3.13, -4.53, -3.41, -4.85, -3.96, -5.49, # 7
					-3.13, -4.68, -3.41, -5.01, -3.96, -5.65, # 8
					-3.13, -4.82, -3.41, -5.15, -3.96, -5.79, # 9
					-3.13, -4.96, -3.41, -5.29, -3.96, -5.94), ncol = 6, byrow = T) # 10
				}
			}
		}
	else if (obs <= 70) {
		# Case I: no intercept, no trend
		if (case == 1) {
			# Doesn't exist! Use asym. crit values
			fnote <- "Small-sample critical values are not provided for Case I. Asymptotic critical values used."
			#		0.10        0.05	   0.010
			# 	 I(0)   I(1)  I(0)  I(1)  I(0) I(1)
			fmat <- matrix(c(3.00, 3.00, 4.20, 4.20, 7.17, 7.17, # 0
				2.44, 3.28, 3.15, 4.11, 4.81, 6.02, # 1
				2.17, 3.19, 2.72, 3.83, 3.88, 5.30, # 2
				2.01, 3.10, 2.45, 3.63, 3.42, 4.84, # 3
				1.90, 3.01, 2.26, 3.48, 3.07, 4.44, # 4
				1.81, 2.93, 2.14, 3.34, 2.82, 4.21, # 5
				1.75, 2.87, 2.04, 3.24, 2.66, 4.05, # 6
				1.70, 2.83, 1.97, 3.18, 2.54, 3.91, # 7
				1.66, 2.79, 1.91, 3.11, 2.45, 3.79, # 8
				1.63, 2.75, 1.86, 3.05, 2.34, 3.68, # 9
				1.60, 2.72, 1.82, 2.99, 2.26, 3.60), ncol = 6, byrow = T) # 10
			if (!is.null(tstat)) {
				tnote <- "Small-sample critical values not provided for Case I. Asymptotic critical values used."
				#		0.10         0.05		   0.010
				# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
				tmat <- matrix(c(-1.62, -1.62, -1.95, -1.95, -2.58, -2.58,# 0
					-1.62, -2.28, -1.95, -2.60, -2.58, -3.22, # 1
					-1.62, -2.68, -1.95, -3.02, -2.58, -3.66, # 2
					-1.62, -3.00, -1.95, -3.33, -2.58, -3.97, # 3
					-1.62, -3.26, -1.95, -3.60, -2.58, -4.23, # 4
					-1.62, -3.49, -1.95, -3.83, -2.58, -4.44, # 5
					-1.62, -3.70, -1.95, -4.04, -2.58, -4.67, # 6
					-1.62, -3.90, -1.95, -4.23, -2.58, -4.88, # 7
					-1.62, -4.09, -1.95, -4.43, -2.58, -5.07, # 8
					-1.62, -4.26, -1.95, -4.61, -2.58, -5.25, # 9
					-1.62, -4.42, -1.95, -4.76, -2.58, -5.44), ncol = 6, byrow = T) # 10
				}
			}
		# Case II: restricted intercept, no trend
		else if (case == 2)	{
			#		0.10         0.05		   0.010
			# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
			fmat <- matrix(c(3.875, 3.875, 4.750, 4.750, 6.740, 6.740, # 0
				3.120, 3.623, 3.780, 4.327, 5.157, 5.957, # 1
				2.730, 3.445, 3.243, 4.043, 4.398, 5.463, # 2
				2.482, 3.310, 2.924, 3.860, 3.916, 5.088, # 3
				2.320, 3.232, 2.725, 3.718, 3.608, 4.860, # 4
				2.193, 3.161, 2.564, 3.650, 3.373, 4.717, # 5
				2.100, 3.121, 2.451, 3.559, 3.180, 4.596, # 6
				2.024, 3.079, 2.351, 3.498, 3.034, 4.426), ncol = 6, byrow = T) # 7
			if (!is.null(tstat)) {
				tnote <- "Critical values do not currently exist for Case II."
				}
			}
		# Case III: unrestricted intercept, no trend
		else if (case == 3)	{
			#		0.10         0.05		   0.010
			# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
			fmat <- matrix(c(6.670, 6.670, 8.370, 8.370, 12.240, 12.240, # 0
				4.125, 4.880, 5.055, 5.915, 7.170, 8.405, # 1
				3.250, 4.237, 3.947, 5.020, 5.487, 6.880, # 2
				2.818, 3.880, 3.370, 4.545, 4.635, 6.055, # 3
				2.552, 3.648, 3.022, 4.256, 4.098, 5.570, # 4
				2.363, 3.510, 2.788, 4.073, 3.747, 5.285, # 5
				2.233, 3.407, 2.629, 3.906, 3.436, 5.044, # 6
				2.138, 3.325, 2.494, 3.786, 3.261, 4.821), ncol = 6, byrow = T) # 7
			if (!is.null(tstat)) {
				tnote <- "Small-sample critical values not provided for Case III. Asymptotic critical values used."
				#		0.10         0.05		   0.010
				# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
				tmat <- matrix(c(-2.57, -2.57, -2.86, -2.86, -3.43, -3.43, # 0
					-2.57, -2.91, -2.86, -3.22, -3.43, -3.82, # 1
					-2.57, -3.21, -2.86, -3.53, -3.43, -4.10, # 2
					-2.57, -3.46, -2.86, -3.78, -3.43, -4.37, # 3
					-2.57, -3.66, -2.86, -3.99, -3.43, -4.60, # 4
					-2.57, -3.86, -2.86, -4.19, -3.43, -4.79, # 5
					-2.57, -4.04, -2.86, -4.38, -3.43, -4.99, # 6
					-2.57, -4.23, -2.86, -4.57, -3.43, -5.19, # 7
					-2.57, -4.40, -2.86, -4.72, -3.43, -5.37, # 8
					-2.57, -4.56, -2.86, -4.88, -3.42, -5.54, # 9
					-2.57, -4.69, 2.86, -5.03, -3.43, -5.68), ncol = 6, byrow = T) # 10
				}
			}
		# Case IV: unrestricted intercept, restricted trend
		else if (case == 4)	{
			#		0.10         0.05		   0.010
			# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
			fmat <- matrix(c(5.530, 5.530, 6.530, 6.530, 8.890, 8.890, # 0
				4.173, 4.647, 4.930, 5.457, 6.577, 7.313, # 1
				3.505, 4.198, 4.100, 4.900, 5.448, 6.435, # 2
				3.098, 3.920, 3.600, 4.512, 4.760, 5.798, # 3
				2.832, 3.738, 3.272, 4.272, 4.293, 5.460, # 4
				2.631, 3.589, 3.043, 4.100, 3.966, 5.234, # 5
				2.485, 3.473, 2.860, 3.951, 3.720, 5.004, # 6
				2.363, 3.394, 2.711, 3.842, 3.509, 4.808), ncol = 6, byrow = T) # 7
			if (!is.null(tstat)) {
				tnote <- "Critical values do not currently exist for Case IV."
				}
			}
		# Case V: unrestricted intercept, unrestricted trend
		else { # case == 5
			#		0.10         0.05		   0.010
			# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
			fmat <- matrix(c(10.020, 10.020, 12.000, 12.000, 16.660, 16.660, # 0
				5.765, 6.455, 6.860, 7.645, 9.370, 10.320, # 1
				4.330, 5.243, 5.110, 6.190, 6.873, 8.163, # 2
				3.615, 4.635, 4.235, 5.363, 5.663, 6.953, # 3
				3.182, 4.258, 3.720, 4.904, 4.922, 6.328, # 4
				2.893, 4.008, 3.368, 4.590, 4.428, 5.898, # 5
				2.683, 3.807, 3.107, 4.343, 4.070, 5.534, # 6
				2.519, 3.669, 2.913, 4.168, 3.774, 5.248), ncol = 6, byrow = T) # 7
			if (!is.null(tstat)) {
				tnote <- "Small-sample critical values not provided for Case V. Asymptotic critical values used."
				#		0.10         0.05		   0.010
				# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
				tmat <- matrix(c(-3.13, -3.13, -3.41, -3.41, -3.96, -3.97, # 0
					-3.13, -3.40, -3.41, -3.69, -3.96, -4.26, # 1
					-3.13, -3.63, -3.41, -3.95, -3.96, -4.53, # 2
					-3.13, -3.84, -3.41, -4.16, -3.96, -4.73, # 3
					-3.13, -4.04, -3.41, -4.36, -3.96, -4.96, # 4
					-3.13, -4.21, -3.41, -4.52, -3.96, -5.13, # 5
					-3.13, -4.37, -3.41, -4.69, -3.96, -5.31, # 6
					-3.13, -4.53, -3.41, -4.85, -3.96, -5.49, # 7
					-3.13, -4.68, -3.41, -5.01, -3.96, -5.65, # 8
					-3.13, -4.82, -3.41, -5.15, -3.96, -5.79, # 9
					-3.13, -4.96, -3.41, -5.29, -3.96, -5.94), ncol = 6, byrow = T) # 10
				}
			}
		}
	else if (obs <= 75) {
		# Case I: no intercept, no trend
		if (case == 1) {
			# Doesn't exist! Use asym. crit values
			fnote <- "Small-sample critical values are not provided for Case I. Asymptotic critical values used."
			#		0.10        0.05	   0.010
			# 	 I(0)   I(1)  I(0)  I(1)  I(0) I(1)
			fmat <- matrix(c(3.00, 3.00, 4.20, 4.20, 7.17, 7.17, # 0
				2.44, 3.28, 3.15, 4.11, 4.81, 6.02, # 1
				2.17, 3.19, 2.72, 3.83, 3.88, 5.30, # 2
				2.01, 3.10, 2.45, 3.63, 3.42, 4.84, # 3
				1.90, 3.01, 2.26, 3.48, 3.07, 4.44, # 4
				1.81, 2.93, 2.14, 3.34, 2.82, 4.21, # 5
				1.75, 2.87, 2.04, 3.24, 2.66, 4.05, # 6
				1.70, 2.83, 1.97, 3.18, 2.54, 3.91, # 7
				1.66, 2.79, 1.91, 3.11, 2.45, 3.79, # 8
				1.63, 2.75, 1.86, 3.05, 2.34, 3.68, # 9
				1.60, 2.72, 1.82, 2.99, 2.26, 3.60), ncol = 6, byrow = T) # 10
			if (!is.null(tstat)) {
				tnote <- "Small-sample critical values not provided for Case I. Asymptotic critical values used."
				#		0.10         0.05		   0.010
				# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
				tmat <- matrix(c(-1.62, -1.62, -1.95, -1.95, -2.58, -2.58, # 0
					-1.62, -2.28, -1.95, -2.60, -2.58, -3.22, # 1
					-1.62, -2.68, -1.95, -3.02, -2.58, -3.66, # 2
					-1.62, -3.00, -1.95, -3.33, -2.58, -3.97, # 3
					-1.62, -3.26, -1.95, -3.60, -2.58, -4.23, # 4
					-1.62, -3.49, -1.95, -3.83, -2.58, -4.44, # 5
					-1.62, -3.70, -1.95, -4.04, -2.58, -4.67, # 6
					-1.62, -3.90, -1.95, -4.23, -2.58, -4.88, # 7
					-1.62, -4.09, -1.95, -4.43, -2.58, -5.07, # 8
					-1.62, -4.26, -1.95, -4.61, -2.58, -5.25, # 9
					-1.62, -4.42, -1.95, -4.76, -2.58, -5.44), ncol = 6, byrow = T) # 10
				}
			}
		# Case II: restricted intercept, no trend
		else if (case == 2)	{
			#		0.10         0.05		   0.010
			# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
			fmat <- matrix(c(3.895, 3.895, 4.760, 4.760, 6.915, 6.915, # 0
				3.133, 3.597, 3.777, 4.320, 5.260, 5.957, # 1
				2.725, 3.455, 3.253, 4.065, 4.458, 5.410, # 2
				2.482, 3.334, 2.946, 3.862, 4.048, 5.092, # 3
				2.313, 3.228, 2.725, 3.718, 3.687, 4.842, # 4
				2.196, 3.166, 2.574, 3.641, 3.427, 4.620, # 5
				2.103, 3.111, 2.449, 3.550, 3.219, 4.526, # 6
				2.023, 3.068, 2.360, 3.478, 3.057, 4.413), ncol = 6, byrow = T) # 7
			if (!is.null(tstat)) {
				tnote <- "Critical values do not currently exist for Case II."
				}
			}
		# Case III: unrestricted intercept, no trend
		else if (case == 3)	{
			#		0.10         0.05		   0.010
			# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
			fmat <- matrix(c(6.720, 6.720, 8.420, 8.420, 12.540, 12.540, # 0
				4.150, 4.885, 5.140, 5.920, 7.225, 8.300, # 1
				3.277, 4.243, 3.983, 5.060, 5.513, 6.860, # 2
				2.838, 3.898, 3.408, 4.550, 4.725, 6.080, # 3
				2.558, 3.654, 3.042, 4.244, 4.168, 5.548, # 4
				2.380, 3.515, 2.802, 4.065, 3.772, 5.213, # 5
				2.244, 3.397, 2.637, 3.900, 3.496, 4.966, # 6
				2.134, 3.313, 2.503, 3.768, 3.266, 4.801), ncol = 6, byrow = T) # 7
			if (!is.null(tstat)) {
				tnote <- "Small-sample critical values not provided for Case III. Asymptotic critical values used."
				#		0.10         0.05		   0.010
				# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
				tmat <- matrix(c(-2.57, -2.57, -2.86, -2.86, -3.43, -3.43, # 0
					-2.57, -2.91, -2.86, -3.22, -3.43, -3.82, # 1
					-2.57, -3.21, -2.86, -3.53, -3.43, -4.10, # 2
					-2.57, -3.46, -2.86, -3.78, -3.43, -4.37, # 3
					-2.57, -3.66, -2.86, -3.99, -3.43, -4.60, # 4
					-2.57, -3.86, -2.86, -4.19, -3.43, -4.79, # 5
					-2.57, -4.04, -2.86, -4.38, -3.43, -4.99, # 6
					-2.57, -4.23, -2.86, -4.57, -3.43, -5.19, # 7
					-2.57, -4.40, -2.86, -4.72, -3.43, -5.37, # 8
					-2.57, -4.56, -2.86, -4.88, -3.42, -5.54, # 9
					-2.57, -4.69, -2.86, -5.03, -3.43, -5.68), ncol = 6, byrow = T) # 10
				}
			}
		# Case IV: unrestricted intercept, restricted trend
		else if (case == 4)	{
			#		0.10         0.05		   0.010
			# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
			fmat <- matrix(c(5.530, 5.530, 6.580, 6.580, 8.905, 8.905, # 0
				4.193, 4.647, 4.937, 5.443, 6.613, 7.253, # 1
				3.505, 4.213, 4.120, 4.855, 5.505, 6.298, # 2
				3.110, 3.900, 3.624, 4.488, 4.808, 5.786, # 3
				2.832, 3.717, 3.298, 4.260, 4.300, 5.377, # 4
				2.636, 3.579, 3.054, 4.079, 3.984, 5.153, # 5
				2.486, 3.469, 2.874, 3.914, 3.728, 4.954, # 6
				2.372, 3.370, 2.718, 3.807, 3.511, 4.789), ncol = 6, byrow = T) # 7
			if (!is.null(tstat)) {
				tnote <- "Critical values do not currently exist for Case IV."
			}
		}
	# Case V: unrestricted intercept, unrestricted trend
	else { # case == 5
			#		0.10         0.05		   0.010
			# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
			fmat <- matrix(c(10.030, 10.030, 12.080, 12.080, 16.610, 16.610, # 0
				5.765, 6.470, 6.880, 7.675, 9.325, 10.325, # 1
				4.323, 5.273, 5.140, 6.153, 6.930, 8.027, # 2
				3.618, 4.630, 4.253, 5.333, 5.698, 6.970, # 3
				3.182, 4.248, 3.724, 4.880, 4.932, 6.224, # 4
				2.890, 3.993, 3.382, 4.567, 4.393, 5.788, # 5
				2.681, 3.800, 3.111, 4.310, 4.060, 5.459, # 6
				2.530, 3.648, 2.915, 4.143, 3.768, 5.229), ncol = 6, byrow = T) # 7
			if (!is.null(tstat)) {
				tnote <- "Small-sample critical values not provided for Case V. Asymptotic critical values used."
				#		0.10         0.05		   0.010
				# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
				tmat <- matrix(c(-3.13, -3.13, -3.41, -3.41, -3.96, -3.97, # 0
					-3.13, -3.40, -3.41, -3.69, -3.96, -4.26, # 1
					-3.13, -3.63, -3.41, -3.95, -3.96, -4.53, # 2
					-3.13, -3.84, -3.41, -4.16, -3.96, -4.73, # 3
					-3.13, -4.04, -3.41, -4.36, -3.96, -4.96, # 4
					-3.13, -4.21, -3.41, -4.52, -3.96, -5.13, # 5
					-3.13, -4.37, -3.41, -4.69, -3.96, -5.31, # 6
					-3.13, -4.53, -3.41, -4.85, -3.96, -5.49, # 7
					-3.13, -4.68, -3.41, -5.01, -3.96, -5.65, # 8
					-3.13, -4.82, -3.41, -5.15, -3.96, -5.79, # 9
					-3.13, -4.96, -3.41, -5.29, -3.96, -5.94), ncol = 6, byrow = T) # 10
				}
			}
		}
	else if (obs <= 80) {
 	# Case I: no intercept, no trend
		if (case == 1) {
			# Doesn't exist! Use asym. crit values
			fnote <- "Small-sample critical values are not provided for Case I. Asymptotic critical values used."
			#		0.10        0.05	   0.010
			# 	 I(0)   I(1)  I(0)  I(1)  I(0) I(1)
			fmat <- matrix(c(3.00, 3.00, 4.20, 4.20, 7.17, 7.17, # 0
				2.44, 3.28, 3.15, 4.11, 4.81, 6.02, # 1
				2.17, 3.19, 2.72, 3.83, 3.88, 5.30, # 2
				2.01, 3.10, 2.45, 3.63, 3.42, 4.84, # 3
				1.90, 3.01, 2.26, 3.48, 3.07, 4.44, # 4
				1.81, 2.93, 2.14, 3.34, 2.82, 4.21, # 5
				1.75, 2.87, 2.04, 3.24, 2.66, 4.05, # 6
				1.70, 2.83, 1.97, 3.18, 2.54, 3.91, # 7
				1.66, 2.79, 1.91, 3.11, 2.45, 3.79, # 8
				1.63, 2.75, 1.86, 3.05, 2.34, 3.68, # 9
				1.60, 2.72, 1.82, 2.99, 2.26, 3.60), ncol = 6, byrow = T) # 10
			if (!is.null(tstat)) {
				tnote <- "Small-sample critical values not provided for Case I. Asymptotic critical values used."
				#		0.10         0.05		   0.010
				# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
				tmat <- matrix(c(-1.62, -1.62, -1.95, -1.95, -2.58, -2.58, # 0
					-1.62, -2.28, -1.95, -2.60, -2.58, -3.22, # 1
					-1.62, -2.68, -1.95, -3.02, -2.58, -3.66, # 2
					-1.62, -3.00, -1.95, -3.33, -2.58, -3.97, # 3
					-1.62, -3.26, -1.95, -3.60, -2.58, -4.23, # 4
					-1.62, -3.49, -1.95, -3.83, -2.58, -4.44, # 5
					-1.62, -3.70, -1.95, -4.04, -2.58, -4.67, # 6
					-1.62, -3.90, -1.95, -4.23, -2.58, -4.88, # 7
					-1.62, -4.09, -1.95, -4.43, -2.58, -5.07, # 8
					-1.62, -4.26, -1.95, -4.61, -2.58, -5.25, # 9
					-1.62, -4.42, -1.95, -4.76, -2.58, -5.44), ncol = 6, byrow = T) # 10
				}
			}
		# Case II: restricted intercept, no trend
		else if (case == 2)	{
			#		0.10         0.05		   0.010
			# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
			fmat <- matrix(c(3.870, 3.870, 4.725, 4.725, 6.695, 6.695, # 0
				3.113, 3.610, 3.740, 4.303, 5.157, 5.917, # 1
				2.713, 3.453, 3.235, 4.053, 4.358, 5.393, # 2
				2.474, 3.312, 2.920, 3.838, 3.908, 5.004, # 3
				2.303, 3.220, 2.688, 3.698, 3.602, 4.787, # 4
				2.303, 3.154, 2.550, 3.606, 3.351, 4.587, # 5
				2.088, 3.103, 2.431, 3.518, 3.173, 4.485, # 6
				2.017, 3.052, 2.336, 3.458, 3.021, 4.350), ncol = 6, byrow = T) # 7
			if (!is.null(tstat)) {
				tnote <- "Critical values do not currently exist for Case II."
				}
			}
		# Case III: unrestricted intercept, no trend
		else if (case == 3)	{
			#		0.10         0.05		   0.010
			# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
			fmat <- matrix(c(6.720, 6.720, 8.400, 8.400, 12.120, 12.120, # 0
				4.135, 4.895, 5.060, 5.930, 7.095, 8.260, # 1
				3.260, 4.247, 3.940, 5.043, 5.407, 6.783, # 2
				2.823, 3.885, 3.363, 4.515, 4.568, 5.960, # 3
				2.548, 3.644, 3.010, 4.216, 4.096, 5.512, # 4
				2.355, 3.500, 2.787, 4.015, 3.725, 5.163, # 5
				2.236, 3.381, 2.627, 3.864, 3.457, 4.943, # 6
				2.129, 3.289, 2.476, 3.746, 3.233, 4.760), ncol = 6, byrow = T) # 7
			if (!is.null(tstat)) {
				tnote <- "Small-sample critical values not provided for Case III. Asymptotic critical values used."
				#		0.10         0.05		   0.010
				# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
				tmat <- matrix(c(-2.57, -2.57, -2.86, -2.86, -3.43, -3.43, # 0
					-2.57, -2.91, -2.86, -3.22, -3.43, -3.82, # 1
					-2.57, -3.21, -2.86, -3.53, -3.43, -4.10, # 2
					-2.57, -3.46, -2.86, -3.78, -3.43, -4.37, # 3
					-2.57, -3.66, -2.86, -3.99, -3.43, -4.60, # 4
					-2.57, -3.86, -2.86, -4.19, -3.43, -4.79, # 5
					-2.57, -4.04, -2.86, -4.38, -3.43, -4.99, # 6
					-2.57, -4.23, -2.86, -4.57, -3.43, -5.19, # 7
					-2.57, -4.40, -2.86, -4.72, -3.43, -5.37, # 8
					-2.57, -4.56, -2.86, -4.88, -3.42, -5.54, # 9
					-2.57, -4.69, -2.86, -5.03, -3.43, -5.68), ncol = 6, byrow = T) # 10
				}
			}
		# Case IV: unrestricted intercept, restricted trend
		else if (case == 4)	{
			#		0.10         0.05		   0.010
			# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
			fmat <- matrix(c(3.870, 3.870, 4.725, 4.725, 6.695, 6.695, # 0
				3.113, 3.610, 3.740, 4.303, 5.157, 5.917, # 1
				2.713, 3.453, 3.235, 4.053, 4.358, 5.393, # 2
				2.474, 3.312, 2.920, 3.838, 3.908, 5.004, # 3
				2.303, 3.220, 2.688, 3.698, 3.602, 4.787, # 4
				2.180, 3.154, 2.550, 3.606, 3.351, 4.587, # 5
				2.088, 3.103, 2.431, 3.518, 3.173, 4.485, # 6
				2.017, 3.052, 2.336, 3.458, 3.021, 4.350), ncol = 6, byrow = T) # 7
			if (!is.null(tstat)) {
				tnote <- "Critical values do not currently exist for Case IV."
				}
			}
		# Case V: unrestricted intercept, unrestricted trend
		else { # case == 5
			#		0.10         0.05		   0.010
			# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
			fmat <- matrix(c(9.960, 9.960, 12.060, 12.060, 16.600, 16.600, # 0
				5.725, 6.450, 6.820, 7.670, 9.170, 10.240, # 1
				4.307, 5.223, 5.067, 6.103, 6.730, 8.053, # 2
				3.588, 4.605, 4.203, 5.320, 5.620, 6.908, # 3
				3.160, 4.230, 3.678, 4.840, 4.890, 6.164, # 4
				2.867, 3.975, 3.335, 4.535, 4.375, 5.703, # 5
				2.657, 3.776, 3.077, 4.284, 4.000, 5.397, # 6
				2.504, 3.631, 2.885, 4.111, 3.728, 5.160), ncol = 6, byrow = T) # 7
			if (!is.null(tstat)) {
				tnote <- "Small-sample critical values not provided for Case V. Asymptotic critical values used."
				#		0.10         0.05		   0.010
				# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
				tmat <- matrix(c(-3.13, -3.13, -3.41, -3.41, -3.96, -3.97, # 0
	 				-3.13, -3.40, -3.41, -3.69, -3.96, -4.26, # 1
					-3.13, -3.63, -3.41, -3.95, -3.96, -4.53, # 2
					-3.13, -3.84, -3.41, -4.16, -3.96, -4.73, # 3
					-3.13, -4.04, -3.41, -4.36, -3.96, -4.96, # 4
					-3.13, -4.21, -3.41, -4.52, -3.96, -5.13, # 5
					-3.13, -4.37, -3.41, -4.69, -3.96, -5.31, # 6
					-3.13, -4.53, -3.41, -4.85, -3.96, -5.49, # 7
					-3.13, -4.68, -3.41, -5.01, -3.96, -5.65, # 8
					-3.13, -4.82, -3.41, -5.15, -3.96, -5.79, # 9
					-3.13, -4.96, -3.41, -5.29, -3.96, -5.94), ncol = 6, byrow = T) # 10
				}
			}
		}
	else { # asymtotic
		# Case I: no intercept, no trend
		if (case == 1) {
			# Doesn't exist! Use asym. crit values
			fnote <- "Asymptotic critical values used."
			#		0.10        0.05	   0.010
			# 	 I(0)   I(1)  I(0)  I(1)  I(0) I(1)
			fmat <- matrix(c(3.00, 3.00, 4.20, 4.20, 7.17, 7.17, # 0
				2.44, 3.28, 3.15, 4.11, 4.81, 6.02, # 1
				2.17, 3.19, 2.72, 3.83, 3.88, 5.30, # 2
				2.01, 3.10, 2.45, 3.63, 3.42, 4.84, # 3
				1.90, 3.01, 2.26, 3.48, 3.07, 4.44, # 4
				1.81, 2.93, 2.14, 3.34, 2.82, 4.21, # 5
				1.75, 2.87, 2.04, 3.24, 2.66, 4.05, # 6
				1.70, 2.83, 1.97, 3.18, 2.54, 3.91, # 7
				1.66, 2.79, 1.91, 3.11, 2.45, 3.79, # 8
				1.63, 2.75, 1.86, 3.05, 2.34, 3.68, # 9
				1.60, 2.72, 1.82, 2.99, 2.26, 3.60), ncol = 6, byrow = T) # 10
			if (!is.null(tstat)) {
				tnote <- "Asymptotic critical values used."
				#		0.10         0.05		   0.010
				# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
				tmat <- matrix(c(-1.62, -1.62, -1.95, -1.95, -2.58, -2.58, # 0
					-1.62, -2.28, -1.95, -2.60, -2.58, -3.22, # 1
					-1.62, -2.68, -1.95, -3.02, -2.58, -3.66, # 2
					-1.62, -3.00, -1.95, -3.33, -2.58, -3.97, # 3
					-1.62, -3.26, -1.95, -3.60, -2.58, -4.23, # 4
					-1.62, -3.49, -1.95, -3.83, -2.58, -4.44, # 5
					-1.62, -3.70, -1.95, -4.04, -2.58, -4.67, # 6
					-1.62, -3.90, -1.95, -4.23, -2.58, -4.88, # 7
					-1.62, -4.09, -1.95, -4.43, -2.58, -5.07, # 8
					-1.62, -4.26, -1.95, -4.61, -2.58, -5.25, # 9
					-1.62, -4.42, -1.95, -4.76, -2.58, -5.44), ncol = 6, byrow = T) # 10
				}
			}
		# Case II: restricted intercept, no trend
		else if (case == 2)	{
			#		0.10         0.05		   0.010
			# 	 I(0)  I(1)  I(0)  I(1)  I(0)  I(1)
			fmat <- matrix(c(3.80, 3.80, 4.60, 4.60, 6.44, 6.44, # 0
				3.02, 3.51, 3.62, 4.16, 4.94, 5.58, # 1
				2.63, 3.35, 3.10, 3.87, 4.13, 5.00, # 2
				2.37, 3.20, 2.79, 3.67, 3.65, 4.66, # 3
				2.20, 3.09, 2.56, 3.49, 3.29, 4.37, # 4
				2.08, 3.00, 2.39, 3.38, 3.06, 4.15, # 5
				1.99, 2.94, 2.27, 3.28, 2.88, 3.99, # 6
				1.92, 2.89, 2.17, 3.21, 2.73, 3.90, # 7
				1.85, 2.85, 2.11, 3.15, 2.62, 3.77, # 8
				1.80, 2.80, 2.04, 3.08, 2.50, 3.68, # 9
				1.76, 2.77, 1.98, 3.04, 2.41, 3.61), ncol = 6, byrow = T) # 10
			fnote <- "Asymptotic critical values used."
			if (!is.null(tstat)) {
				tnote <- "Critical values do not currently exist for Case II."
				}
			}
		# Case III: unrestricted intercept, no trend
		else if (case == 3)	{
			#		0.10         0.05		   0.010
			# 	 I(0)  I(1)  I(0)  I(1)  I(0)  I(1)
			fmat <- matrix(c(6.58, 6.58, 8.21, 8.21, 11.79, 11.79, # 0
				4.04, 4.78, 4.94, 5.73, 6.84, 7.84, # 1
				3.17, 4.14, 3.79, 4.85, 5.15, 6.36, # 2
				2.72, 3.77, 3.23, 4.35, 4.29, 5.61, # 3
				2.45, 3.52, 2.86, 4.01, 3.74, 5.06, # 4
				2.26, 3.35, 2.62, 3.79, 3.41, 4.68, # 5
				2.12, 3.23, 2.45, 3.61, 3.15, 4.43, # 6
				2.03, 3.13, 2.32, 3.50, 2.96, 4.26, # 7
				1.95, 3.06, 2.22, 3.39, 2.79, 4.10, # 8
				1.88, 2.99, 2.14, 3.30, 2.65, 3.97, # 9
				1.83, 2.94, 2.06, 3.24, 2.54, 3.86), ncol = 6, byrow = T) # 10
			fnote <- "Asymptotic critical values used."
			if (!is.null(tstat)) {
				tnote <- "Asymptotic critical values used."
				#		0.10         0.05		   0.010
				# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
				tmat <- matrix(c(-2.57, -2.57, -2.86, -2.86, -3.43, -3.43, # 0
					-2.57, -2.91, -2.86, -3.22, -3.43, -3.82, # 1
					-2.57, -3.21, -2.86, -3.53, -3.43, -4.10, # 2
					-2.57, -3.46, -2.86, -3.78, -3.43, -4.37, # 3
					-2.57, -3.66, -2.86, -3.99, -3.43, -4.60, # 4
					-2.57, -3.86, -2.86, -4.19, -3.43, -4.79, # 5
					-2.57, -4.04, -2.86, -4.38, -3.43, -4.99, # 6
					-2.57, -4.23, -2.86, -4.57, -3.43, -5.19, # 7
					-2.57, -4.40, -2.86, -4.72, -3.43, -5.37, # 8
					-2.57, -4.56, -2.86, -4.88, -3.42, -5.54, # 9
					-2.57, -4.69, -2.86, -5.03, -3.43, -5.68), ncol = 6, byrow = T) # 10
				}
			}
		# Case IV: unrestricted intercept, restricted trend
		else if (case == 4)	{
			#		0.10         0.05		   0.010
			# 	 I(0)  I(1)  I(0)  I(1)  I(0)  I(1)
			fmat <- matrix(c(5.37, 5.37, 6.29, 6.29, 8.26, 8.26, # 0
				4.05, 4.49, 4.68, 5.15, 6.10, 6.73, # 1
				3.38, 4.02, 3.88, 4.61, 4.99, 5.85, # 2
				2.97, 3.74, 3.38, 4.23, 4.30, 5.23, # 3
				2.68, 3.53, 3.05, 3.97, 3.81, 4.92, # 4
				2.49, 3.38, 2.81, 3.76, 3.50, 4.63, # 5
				2.33, 3.25, 2.63, 3.62, 3.27, 4.39, # 6
				2.22, 3.17, 2.50, 3.50, 3.07, 4.23, # 7
				2.13, 3.09, 2.38, 3.41, 2.93, 4.06, # 8
				2.05, 3.02, 2.30, 3.33, 2.79, 3.93, # 9
				1.98, 2.97, 2.21, 3.25, 2.68, 3.84), ncol = 6, byrow = T) # 10
			fnote <- "Asymptotic critical values used."
			if (!is.null(tstat)) {
				tnote <- "Critical values do not currently exist for Case IV."
				}
			}
		# Case V: unrestricted intercept, unrestricted trend
		else { # case == 5
			#		0.10         0.05		   0.010
			# 	 I(0)  I(1)  I(0)  I(1)  I(0)  I(1)
			fmat <- matrix(c(9.81, 9.81, 11.64, 11.64, 15.73, 15.73, # 0
				5.59, 6.26, 6.56, 7.30, 8.74, 9.63, # 1
				4.19, 5.06, 4.87, 5.85, 6.34, 7.52, # 2
				3.47, 4.45, 4.01, 5.07, 5.17, 6.36, # 3
				3.03, 4.06, 3.47, 4.57, 4.40, 5.72, # 4
				2.75, 3.79, 3.12, 4.25, 3.93, 5.23, # 5
				2.53, 3.59, 2.87, 4.00, 3.60, 4.90, # 6
				2.38, 3.45, 2.69, 3.83, 3.34, 4.63, # 7
				2.26, 3.34, 2.55, 3.68, 3.15, 4.43, # 8
				2.16, 3.24, 2.43, 3.56, 2.97, 4.24, # 9
				2.07, 3.16, 2.33, 3.46, 2.84, 4.10), ncol = 6, byrow = T) # 10
			fnote <- "Asymptotic critical values used."
			if (!is.null(tstat)) {
				tnote <- "Asymptotic critical values used."
				#		0.10         0.05		   0.010
				# 	 I(0)   I(1)    I(0)   I(1)  I(0)    I(1)
				tmat <- matrix(c(-3.13, -3.13, -3.41, -3.41, -3.96, -3.97, # 0
					-3.13, -3.40, -3.41, -3.69, -3.96, -4.26, # 1
					-3.13, -3.63, -3.41, -3.95, -3.96, -4.53, # 2
					-3.13, -3.84, -3.41, -4.16, -3.96, -4.73, # 3
					-3.13, -4.04, -3.41, -4.36, -3.96, -4.96, # 4
					-3.13, -4.21, -3.41, -4.52, -3.96, -5.13, # 5
					-3.13, -4.37, -3.41, -4.69, -3.96, -5.31, # 6
					-3.13, -4.53, -3.41, -4.85, -3.96, -5.49, # 7
					-3.13, -4.68, -3.41, -5.01, -3.96, -5.65, # 8
					-3.13, -4.82, -3.41, -5.15, -3.96, -5.79, # 9
					-3.13, -4.96, -3.41, -5.29, -3.96, -5.94), ncol = 6, byrow = T) # 10
				}
			}
		}
	# We have an fmat and a tmat from above, as well as potentially a fnote and tnote
	# Grabbing crit values from fmat/tmat above
	# Narayan only provides values for k <= 7. If 7 < k, we're using k = 7 but noting it.
	# For asymptotic: if k > 10, we're using k = 10 and noting it.
	f_10_1 <- f_10_1 <- f_05_1 <- f_05_1 <- f_01_1 <- f_01_1 <- NULL
	if (obs <= 80) { # max rows k = 7
		if (k < 7) {
			k2 <- k + 1 # b/c of 0'th row
			f_10_0 <- fmat[k2,1]
			f_10_1 <- fmat[k2,2]
			f_05_0 <- fmat[k2,3]
			f_05_1 <- fmat[k2,4]
			f_01_0 <- fmat[k2,5]
			f_01_1 <- fmat[k2,6]
			}
		else {
			f_10_0 <- fmat[8,1]
			f_10_1 <- fmat[8,2]
			f_05_0 <- fmat[8,3]
			f_05_1 <- fmat[8,4]
			f_01_0 <- fmat[8,5]
			f_01_1 <- fmat[8,6]
			if (k > 7) {
				fnote <- paste(fnote, "Small-sample critical values only available up to k = 7.", sep = " ")
				}
			}
		}
	else { # asym max rows k=10
		if (k < 10) {
			k2 <- k + 1
			f_10_0 <- fmat[k2,1]
			f_10_1 <- fmat[k2,2]
			f_05_0 <- fmat[k2,3]
			f_05_1 <- fmat[k2,4]
			f_01_0 <- fmat[k2,5]
			f_01_1 <- fmat[k2,6]
			}
		else {
			f_10_0 <- fmat[11,1]
			f_10_1 <- fmat[11,2]
			f_05_0 <- fmat[11,3]
			f_05_1 <- fmat[11,4]
			f_01_0 <- fmat[11,5]
			f_01_1 <- fmat[11,6]
			if (k > 10)	{
				fnote <- paste(fnote, "Asymptotic critical values only available up to k = 10.")
				}
			}
		}
	# Initialize some output
	toutput <- ""
	if (!is.null(tstat)) {
		if (case == 2 | case == 4) {
			}
		else {
			if (k < 10)	{
				k2 <- k + 1
				t_10_0 <- tmat[k2,1]
				t_10_1 <- tmat[k2,2]
				t_05_0 <- tmat[k2,3]
				t_05_1 <- tmat[k2,4]
				t_01_0 <- tmat[k2,5]
				t_01_1 <- tmat[k2,6]
				}
			else {
				t_10_0 <- tmat[11,1]
				t_10_1 <- tmat[11,2]
				t_05_0 <- tmat[11,3]
				t_05_1 <- tmat[11,4]
				t_01_0 <- tmat[11,5]
				t_01_1 <- tmat[11,6]
				tnote <- paste(tnote, "Asymptotic critical values only available up to k = 10.")
				}
			toutput <- paste("------------------------------------------------------",
				"-                       t-test                       -",
				"------------------------------------------------------",
				"                <------- I(0) ------------ I(1) ----->",
				paste("10% critical value", paste( format(t_10_0, digits = 3, nsmall = 2), format(t_10_1, digits = 3, nsmall = 2), sep = "            "), sep = "       "),
				paste("5% critical value", paste( format(t_05_0, digits = 3, nsmall = 2), format(t_05_1, digits = 3, nsmall = 2), sep = "            "), sep = "        "),
				paste("1% critical value", paste( format(t_01_0, digits = 3, nsmall = 2), format(t_01_1, digits = 3, nsmall = 2), sep = "            "), sep = "        "),
				"\n",
				paste("t statistic = ", tstat, sep = ""), sep = "\n ")
			}
		}
	# Output
	if(case == 1) {
		cnote <- "(No intercept; no trend)"
	} else if(case == 2) {
		cnote <- "(Intercept included in F-stat restriction; no trend)"		
	} else if(case == 3) {
		cnote <- "(Unrestricted intercept; no trend)"		
	} else if(case == 4) {
		cnote <- "(Unrestricted intercept; trend included in F-stat restriction)"		
	} else if(case == 5) {
		cnote <- "(Unrestricted intercept; unrestricted trend)"		
	}
	cat("\n",
		"PESARAN, SHIN AND SMITH (2001) COINTEGRATION TEST", "\n\n",
		paste("Observations: ", obs, sep = ""), "\n",
		paste("Number of Lagged Regressors (not including LDV) (k): ", k, sep = ""), "\n",
		paste("Case:", case, cnote, sep = " "), "\n\n",
		"------------------------------------------------------", "\n",
		"-                       F-test                       -", "\n",
		"------------------------------------------------------", "\n",
		"                <------- I(0) ------------ I(1) ----->", "\n",
		paste("10% critical value", paste(format(f_10_0, digits = 3, nsmall = 2), format(f_10_1, digits = 3, nsmall = 2), sep = "            "), sep = "       "), "\n",
		paste("5% critical value", paste(format(f_05_0, digits = 3, nsmall = 2), format(f_05_1, digits = 3, nsmall = 2), sep = "            "), sep = "        "), "\n",
		paste("1% critical value", paste(format(f_01_0, digits = 3, nsmall = 2), format(f_01_1, digits = 3, nsmall = 2), sep = "            "), sep = "        "), "\n",
		"\n\n",
		paste("F-statistic = ", fstat, sep = ""), "\n",
		toutput,
		"\n",
		"------------------------------------------------------", "\n",
		if(!is.null(fnote)) {paste("F-statistic note: ", fnote, sep = "")}, "\n",
		if(!is.null(tnote)) {paste("t-statistic note: ", tnote, sep = "")}, "\n")
	if(toutput != "") {
		out <- structure(list(k = k, obs = obs, fstat = fstat, tstat = tstat, case = case,
			ftest.I0.p10 = f_10_0, ftest.I1.p10 = f_10_1,
			ftest.I0.p05 = f_05_0, ftest.I1.p05 = f_05_1,
			ftest.I0.p01 = f_01_0, ftest.I1.p01 = f_01_1,
			ttest.I1.p10 = t_10_1, ttest.I0.p10 = t_10_0,
			ttest.I1.p05 = t_05_1, ttest.I0.p05 = t_05_0,
			ttest.I1.p01 = t_01_1, ttest.I0.p01 = t_01_0))
		}
	else {
		out <- structure(list(k = k, obs = obs, fstat = fstat, tstat = tstat, case = case,
			ftest.I0.p10 = f_10_0, ftest.I1.p10 = f_10_1,
			ftest.I0.p05 = f_05_0, ftest.I1.p05 = f_05_1,
			ftest.I0.p01 = f_01_0, ftest.I1.p01 = f_01_1))
	}
	if(object.out == TRUE) {
		out
	}	
}


##########################################
# -----(8) dynardl.auto.correlated ------#
##########################################
#' Run a variety of autocorrelation tests on the residuals from a \code{\link{dynardl}} model
#' @param x a \code{dynardl} model
#' @param bg.type a character string for the type of Breusch-Godfrey test to run. The default is \code{Chisq}: the Chisq test statistic. The other option is \code{F}: the F-test statistic
#' @param digits the number of digits to round to when showing output. The default is \code{3}
#' @param order the maximum order of serial autocorrelation to test when executing the Breusch-Godfrey test
#' @param object.out if \code{TRUE}, and \code{dynardl.auto.correlated} is assigned to an object, the AIC, BIC, and results will be stored for the user's convenience
#' @return The results of autocorrelation tests
#' @details
#' This is a simple and convenient way to test whether the residuals from the \code{dynardl} model are white noise. As an aside, this is also why \code{dynardl} has a \code{simulate = FALSE} argument: users can ensure the model has white noise residuals before estimating a potentially time-intensive simulation. The output also reminds the user of the null hypotheses for the autocorrelation tests
#' @importFrom lmtest bgtest
#' @importFrom stats shapiro.test
#' @author Soren Jordan and Andrew Q. Philips
#' @keywords utilities
#' @examples
#' # Using the ineq data from dynamac
#' ardl.model <- dynardl(concern ~ incshare10 + urate, data = ineq, 
#'        lags = list("concern" = 1, "incshare10" = 1),
#'        diffs = c("incshare10", "urate"), 
#'        lagdiffs = list("concern" = 1),
#'        ec = TRUE, simulate = FALSE)
#' dynardl.auto.correlated(ardl.model)
#' @export

dynardl.auto.correlated <- function(x, bg.type = "Chisq", digits = 3, order = NULL, object.out = FALSE) { # dw.alt = "greater",
	out <- list()
	if(!(identical(class(x), "dynardl"))) { # Check if data are from dynardl
		stop("To test autocorrelation, provide a dynardl model object.")
	}
	if(!(is.null(order))) {
		out$bg <- bgtest(x$model, type = bg.type, order = order)
	} else {
		out$bg <- bgtest(x$model, type = bg.type)
	}
	# out$dw <- dwtest(x$model, alternative = dw.alt)
	out$sw <- shapiro.test(x$model$residuals)
	out$AIC <- round(AIC(x$model), digits = digits)
	out$BIC <- round(BIC(x$model), digits = digits)
	out$logLik <- round(logLik(x$model), digits = digits)
	flush.console()
	cat("\n",
		"------------------------------------------------------", "\n",
		"Breusch-Godfrey LM Test", "\n",
		paste("Test statistic:", round(out$bg$statistic, digits = digits), sep = " "), "\n",
		paste("p-value:", round(out$bg$p.value, digits = digits), sep = " "), "\n",
		paste("H_0: no autocorrelation up to AR", out$bg$parameter[1], sep = " "), "\n",
		
		"\n",
		"------------------------------------------------------", "\n",
		#
		#"Durbin-Watson Test", "\n",
		#paste("Test statistic:", round(out$dw$statistic, digits = 3), sep = " "), "\n",
		#paste("p-value:", round(out$dw$p.value, digits = 3), sep = " "), "\n",
		#paste("H_0: no AR(1) autocorrelation"), "\n",
		#if(dw.alt == "greater") {
		#print("Durbin-Watson null hypothesis: autocorrelation of residuals is 0. Alternative hypothesis: autocorrelation of residuals is greater than 0.")	
		#} else if(dw.alt == "two.sided") {
		#	print("Durbin-Watson null hypothesis: autocorrelation of residuals is 0. Alternative hypothesis: autocorrelation of residuals is greater OR less than 0 (two-sided).")		
		#} else if(dw.alt == "less") {
		#print("Durbin-Watson null hypothesis: autocorrelation of residuals is 0. Alternative hypothesis: autocorrelation of residuals is less than 0.")		
		#}
		#	
		#"\n",
		#"------------------------------------------------------", "\n",
		
		"Shapiro-Wilk Test for Normality", "\n",
		paste("Test statistic:", round(out$sw$statistic, digits = digits), sep = " "), "\n",
		paste("p-value:", round(out$sw$p.value, digits = digits), sep = " "), "\n",
		paste("H_0: residuals are distributed normal"), "\n",
		
		"\n",
		"------------------------------------------------------", "\n",
		paste("Log-likelihood:", round(out$logLik, digits = digits), sep = " "), "\n",
		paste("AIC:", round(out$AIC, digits = digits), sep = " "), "\n",
		paste("BIC:", round(out$BIC, digits = digits), sep = " "), "\n",
		paste("Note: AIC and BIC calculated with k =", x$model$rank, "on T =", length(x$model$residuals), "observations.", sep = " "), "\n",
		
		"\n",
		"------------------------------------------------------", "\n")
		if(round(out$bg$p.value, digits = digits) < 0.01) {
			cat("Breusch-Godfrey test indicates we reject the null hypothesis of no autocorrelation at p < 0.01.", "\n",	
			"Add lags to remove autocorrelation before running dynardl simulations.", "\n")	
		} else if(round(out$bg$p.value, digits = digits) < 0.05) {
			cat("Breusch-Godfrey test indicates we reject the null hypothesis of no autocorrelation at p < 0.05.", "\n",
			"Add lags to remove autocorrelation before running dynardl simulations.", "\n")	
		} else if(round(out$bg$p.value, digits = digits) < 0.10) {
			cat("Breusch-Godfrey test indicates we reject the null hypothesis of no autocorrelation at p < 0.10.", "\n", 		
			"Add lags to remove autocorrelation before running dynardl simulations.", "\n")	
		}
		if(round(out$sw$p.value, digits = digits) < 0.01) {
			cat("Shapiro-Wilk test indicates we reject the null hypothesis of normality at p < 0.01.", "\n")	
		} else if(round(out$sw$p.value, digits = digits) < 0.05) {
			cat("Shapiro-Wilk test indicates we reject the null hypothesis of normality at p < 0.05.", "\n")		
		} else if(round(out$sw$p.value, digits = digits) < 0.10) {
			cat("Shapiro-Wilk test indicates we reject the null hypothesis of normality at p < 0.10.", "\n")		
		}
		if(object.out == TRUE) {
			out
	}		
}



##################################
# -----(9) summary.dynardl ------#
##################################
#' Enable summary calls to \code{\link{dynardl}} model objects
#' @param object a \code{dynardl} model
#' @param ... additional arguments in the generic summary call
#' @return A summary of the fitted ARDL model.
#' @details
#' \code{dynardl}, by default, stores regression results in \code{foo$model}. This calls those results directly with \code{summary}
#' @author Soren Jordan and Andrew Q. Philips
#' @keywords utilities
#' @examples
#' # Using the ineq data from dynamac
#' ardl.model <- dynardl(concern ~ incshare10 + urate, data = ineq, 
#'        lags = list("concern" = 1, "incshare10" = 1),
#'        diffs = c("incshare10", "urate"), 
#'        lagdiffs = list("concern" = 1),
#'        ec = TRUE, simulate = FALSE)
#' summary(ardl.model)
#' @export

summary.dynardl <- function(object, ...) {
	stopifnot(inherits(object, "dynardl"))
	summary(object$model)
}

#############################################
# ------(10) dynardl.simulation.plot -------#
#############################################
#' Create a plot of a simulated response in a \code{\link{dynardl}} model
#' @param x a \code{dynardl} model with a simulation to be plotted
#' @param type whether the plot should be an area plot (\code{area}) or a spike plot (\code{spike})
#' @param response whether the plot of the response should be shown in levels of the dependent variable (\code{levels}), levels from the mean of the dependent variable (\code{levels.from.mean}), period-over-period changes in the dependent variable (\code{diffs}), the absolute value of the (decreasing) change in the dependent variable  in each time period due to the shock (\code{shock.effect.decay}), the sum of the period-over-period changes (\code{cumulative.diffs}), or the absolute value of the cumulative differences (where negative effects are treated as positive) (\code{cumulative.abs.diffs}). The default is \code{levels}
#' @param bw should the colors be in black and white (for publication)? The default is \code{FALSE}
#' @param ylab a user-defined y-label to be used instead of the default
#' @param xlab a user-defined x-label to be used instead of the default
#' @param ylim a user-defined y-limit to be used instead of the default (for instance, for shared axes)
#' @param tol when deciding when to stop calculating the absolute value of the shocks to the dependent variable, you can specify the minimum amount of movement required to qualify as a non-noise change over time periods (for calculating absolute cumulative differences). The default is 0.1 percent of the mean of the dependent variable. Specify a \code{tol} or a \code{last.period}. If both are specified, \code{last.period} overrides \code{tol}
#' @param start.period which period of the simulation to begin the plot with. You can view the equilibriating behavior of the dependent variable, or you can skip forward in time (maybe to just before the shock). The default is \code{1} (the first period of the simulation)
#' @param last.period when deciding when to stop calculating the absolute value of the shocks to the dependent variable, you can specify a specific period in which to stop calculating absolute cumulative differences. Specify a \code{tol} or a \code{last.period}. If both are specified, \code{last.period} overrides \code{tol}
#' @param abs.errors when calculating confidence for the absolute cumulative effect, should differences accumulate in each time time period (\code{cumulate}, which could be explosive if the error in the model is large), should differences be observed at each time (\code{within.period}, which will have smaller values in equilibrium than when changing), or should only the values be plotted (\code{none}). The default is \code{none}
#' @param ... other arguments to be passed to the call to plot
#' @return a plot of the simulated dynardl model
#' @details
#' When running \code{dynardl}, \code{simulate} must be \code{TRUE} so that there is a simulation to plot. For types \code{cumulative.diffs} and \code{cumulative.abs.diffs}, \code{fullsims} must be \code{TRUE} in the \code{dynardl} simulation
#' @importFrom graphics lines plot points polygon segments
#' @author Soren Jordan and Andrew Q. Philips
#' @keywords utilities
#' @examples
#' # Using the ineq data in dynamac
#' # Shocking Income Top 10
#' set.seed(1)
#' ardl.model <- dynardl(concern ~ incshare10 + urate, data = ineq, 
#'        lags = list("concern" = 1, "incshare10" = 1),
#'        diffs = c("incshare10", "urate"), 
#'        lagdiffs = list("concern" = 1),
#'        ec = TRUE, simulate = TRUE, range = 30,
#'        shockvar = "incshare10", fullsims = TRUE)
#' 
#' # Shows absolute levels
#' dynardl.simulation.plot(ardl.model)	
#' # Shows changes from mean level
#' dynardl.simulation.plot(ardl.model, response = "levels.from.mean")  
#' # Same plot, but with spikeplot
#' dynardl.simulation.plot(ardl.model, type = "spike", response = "levels.from.mean")  
#' # Grayscale plots
#' dynardl.simulation.plot(ardl.model, bw = TRUE)	 
#' @export

dynardl.simulation.plot <- function(x, type = "area", response = "levels", bw = FALSE, last.period = NULL, tol = (abs(x$model$ymean) * 0.01), start.period = 1, abs.errors = "none", ylim = NULL, ylab = NULL, xlab = NULL, ...) {
	if(!(identical(class(x), "dynardl"))) { # Check if data are from dynardl
		stop("To plot simulation, provide a dynardl model object (with a simulation).")
	}
	if(x$model$simulate == FALSE) {
		stop("dynardl object does not include simulation to plot.")
	}
	if(!(response %in% c("levels", "levels.from.mean", "diffs", "cumulative.diffs", "cumulative.abs.diffs", "shock.effect.decay"))) {
		stop("Response must be one of 'levels', 'levels.from.mean', 'diffs', 'shock.effect.decay', 'cumulative.diffs', or 'cumulative.abs.diffs'.")
	}
	if(!(type %in% c("area", "spike"))) {
		stop("Plot type must be either an area plot ('area') or a spike plot ('spike').")
	}
	if(x$model$ldv == FALSE) {
		warning("Responses are from a model with no lagged dependent variable (LDV): are you sure you want this?")
	}
	z <- x$simulation
	z <- z[c(start.period:length(z$shocktime)),]
	if(!(identical(last.period, NULL))) {
		if(last.period > length(x$simulation$central)) {
			warning(paste(paste("last.period requested exceeds simulation range. Calculating on"), paste(length(x$simulation$central)), paste("periods."), sep = " "))
			is.changing <- last.period <- length(x$simulation$central)
		}
		warning(paste(paste("Movement in Y tested on "), paste(last.period), paste(" periods: regardless of if Y might still be moving OR plot might include noise in cumulative absolute difference."), sep = ""))
		is.changing <- last.period
	} else { # If we're calculating based on tolerance
		is.changing <- NULL
		for(i in 1:length(x$simulation$d.central)) { # Test if it's changing
			if(abs(x$simulation$d.central[i]) > tol) {
				is.changing <- i
			}
		}
	}
	if(response == "levels") { # If we're just plotting levels of Y
		plot(z$time, z$ll95, type = "n", 
			ylim = ifelse(c(is.null(ylim), is.null(ylim)), c(min(z$ll95), max(z$ul95)), ylim),
			ylab = ifelse(is.null(ylab), "Y Value", ylab), xlab = ifelse(is.null(xlab), "Time", xlab), ...)
		if(type == "area") { 
			polygon(c(z$time, rev(z$time)), c(z$ul95, rev(z$ll95)), col = ifelse(bw == FALSE, "skyblue1", "grey70"), border = NA) # 95
			polygon(c(z$time, rev(z$time)), c(z$ul90, rev(z$ll90)), col = ifelse(bw == FALSE, "skyblue3", "grey50"), border = NA) # 90
			polygon(c(z$time, rev(z$time)), c(z$ul75, rev(z$ll75)), col = "grey30", border = NA) # 75
			# Actual response
			lines(z$time, z$central, lty = 2, lwd = 3)
		} else { # if it's a spikeplot
			for(i in 1:length(z$time)) { # 95 percent sig
				segments(z$time[i], z$ll95[i], z$time[i], z$ul95[i], lwd = 1, col = ifelse(bw == FALSE, "skyblue1", "grey70"))
			}
			for(i in 1:length(z$time)) { # 90 percent sig
				segments(z$time[i], z$ll90[i], z$time[i], z$ul90[i], lwd = 3, col = ifelse(bw == FALSE, "skyblue3", "grey50"))
			}
			for(i in 1:length(z$time)) { # 75 percent sig
				segments(z$time[i], z$ll75[i], z$time[i], z$ul75[i], lwd = 5, col = "grey30")
			}
			# Actual response
			points(z$time, z$central, lwd = 4)	
		}
	} 
	else if(response == "levels.from.mean") { # If it's changes from the mean, changes values, same code
		z$ll95 <- z$ll95 - x$model$ymean
		z$ll90 <- z$ll90 - x$model$ymean
		z$ll75 <- z$ll75 - x$model$ymean
		z$ul75 <- z$ul75 - x$model$ymean
		z$ul90 <- z$ul90 - x$model$ymean
		z$ul95 <- z$ul95 - x$model$ymean
		z$central <- z$central - x$model$ymean
		plot(z$time, z$ll95, type = "n",
			ylim = ifelse(c(is.null(ylim), is.null(ylim)), c(min(z$ll95), max(z$ul95)), ylim),
			ylab = ifelse(is.null(ylab), "Changes from Y Mean Value", ylab), xlab = ifelse(is.null(xlab), "Time", xlab), ...)
		if(type == "area") {
			polygon(c(z$time, rev(z$time)), c(z$ul95, rev(z$ll95)), col = ifelse(bw == FALSE, "skyblue1", "grey70"), border = NA) # 95
			polygon(c(z$time, rev(z$time)), c(z$ul90, rev(z$ll90)), col = ifelse(bw == FALSE, "skyblue3", "grey50"), border = NA) # 90
			polygon(c(z$time, rev(z$time)), c(z$ul75, rev(z$ll75)), col = "grey30", border = NA) # 75
			# Actual response
			lines(z$time, z$central, lty = 2, lwd = 3)
		} else { # if it's a spikeplot
			for(i in 1:length(z$time)) { # 95 percent sig
				segments(z$time[i], z$ll95[i], z$time[i], z$ul95[i], lwd = 1, col = ifelse(bw == FALSE, "skyblue1", "grey70"))
			}
			for(i in 1:length(z$time)) { # 90 percent sig
				segments(z$time[i], z$ll90[i], z$time[i], z$ul90[i], lwd = 3, col = ifelse(bw == FALSE, "skyblue3", "grey50"))
			}
			for(i in 1:length(z$time)) { # 75 percent sig
				segments(z$time[i], z$ll75[i], z$time[i], z$ul75[i], lwd = 5, col = "grey30")
			}
			# Actual response
			points(z$time, z$central, lwd = 4)
		}
	}
	else if(response == "diffs") { # If it's differences in Y
		plot(z$time, z$d.ll95, type = "n",
			ylim = ifelse(c(is.null(ylim), is.null(ylim)), c(min(z$d.ll95), max(z$d.ul95)), ylim),
			ylab = ifelse(is.null(ylab), "Change in Y Value", ylab), xlab = ifelse(is.null(xlab), "Time", xlab), ...)
		if(type == "area") { 
			polygon(c(z$time, rev(z$time)), c(z$d.ul95, rev(z$d.ll95)), col = ifelse(bw == FALSE, "skyblue1", "grey70"), border = NA) # 95
			polygon(c(z$time, rev(z$time)), c(z$d.ul90, rev(z$d.ll90)), col = ifelse(bw == FALSE, "skyblue3", "grey50"), border = NA) # 90
			polygon(c(z$time, rev(z$time)), c(z$d.ul75, rev(z$d.ll75)), col = "grey30", border = NA) # 75
			# Actual response
			lines(z$time, z$d.central, lty = 2, lwd = 3)
		} else { # if it's a spikeplot
			for(i in 1:length(z$time)) { # 95 percent sig
				segments(z$time[i], z$d.ll95[i], z$time[i], z$d.ul95[i], lwd = 1, col = ifelse(bw == FALSE, "skyblue1", "grey70"))
			}
			for(i in 1:length(z$time)) { # 90 percent sig
				segments(z$time[i], z$d.ll90[i], z$time[i], z$d.ul90[i], lwd = 3, col = ifelse(bw == FALSE, "skyblue3", "grey50"))
			}
			for(i in 1:length(z$time)) { # 75 percent sig
				segments(z$time[i], z$d.ll75[i], z$time[i], z$d.ul75[i], lwd = 5, col = "grey30")
			}
			# Actual response
			points(z$time, z$d.central, lwd = 4)	
		}
	}
	else if(response == "cumulative.diffs") { # If it's cumulative differences in Y
		if(identical(x$rawsims, NULL)) {
			stop("dynardl object must have fullsims = TRUE to track cumulative differences.")
		}
		a <- x$rawsims[,start.period:ncol(x$rawsims)]
		# Two matrices: one for the diffs of the raw simulations (by simulation), the other for cumulative
		# Here, diffs sims is going to be 2:length of simulation, with the first period NA. This is to sync it with the d.central output
		diff.sims <- cum.diff.sims <- matrix(rep(NA, nrow(a)*(ncol(a) - 1)), nrow = nrow(a)) # Last column is central tendency
		# Track the differences for each simulation
		for(i in 2:ncol(diff.sims)) {
			diff.sims[,i] <- a[,i] - a[,(i - 1)]
		}
		# First period: no cumulative difference
		cum.diff.sims[,2] <- diff.sims[,2]
		for(i in 3:ncol(cum.diff.sims)) {
			cum.diff.sims[,i] <- rowSums(diff.sims[,1:i], na.rm = T)
		}
		d.ll95 <- d.ul95 <- d.ll90 <- d.ul90 <- d.ll75 <- d.ul75 <- d.central <- rep(NA, ncol(cum.diff.sims))
		for(i in 1:ncol(cum.diff.sims)) {
			d.ll95[i] <- quantile(cum.diff.sims[,i], 0.025, na.rm = T)
			d.ll90[i] <- quantile(cum.diff.sims[,i], 0.050, na.rm = T)
			d.ll75[i] <- quantile(cum.diff.sims[,i], 0.125, na.rm = T)
			d.ul75[i] <- quantile(cum.diff.sims[,i], 0.875, na.rm = T)
			d.ul90[i] <- quantile(cum.diff.sims[,i], 0.950, na.rm = T)
			d.ul95[i] <- quantile(cum.diff.sims[,i], 0.975, na.rm = T)
			d.central[i] <- ifelse(identical(x$rawsims$central[1], "median"), median(cum.diff.sims[,i], na.rm = T), mean(cum.diff.sims[,i], na.rm = T))
		}
		time <- seq(start.period, (ncol(cum.diff.sims) + start.period - 1), 1)
		plot(time, d.ll95, type = "n", 
			ylim = ifelse(c(is.null(ylim), is.null(ylim)), c(min(d.ll95, na.rm = T), max(d.ul95, na.rm = T)), ylim),
			ylab = ifelse(is.null(ylab), "Cumulative Change in Y Value", ylab), xlab = ifelse(is.null(xlab), "Time", xlab), ...)
		if(type == "area") { 
			polygon(c(time, rev(time)), c(d.ul95, rev(d.ll95)), col = ifelse(bw == FALSE, "skyblue1", "grey70"), border = NA) # 95
			polygon(c(time, rev(time)), c(d.ul90, rev(d.ll90)), col = ifelse(bw == FALSE, "skyblue3", "grey50"), border = NA) # 90
			polygon(c(time, rev(time)), c(d.ul75, rev(d.ll75)), col = "grey30", border = NA) # 75
			# Actual response
			lines(time, d.central, lty = 2, lwd = 3)
		} else { # if it's a spikeplot
			for(i in 1:length(time)) { # 95 percent sig
				segments(time[i], d.ll95[i], time[i], d.ul95[i], lwd = 1, col = ifelse(bw == FALSE, "skyblue1", "grey70"))
			}
			for(i in 1:length(time)) { # 90 percent sig
				segments(time[i], d.ll90[i], time[i], d.ul90[i], lwd = 3, col = ifelse(bw == FALSE, "skyblue3", "grey50"))
			}
			for(i in 1:length(z$time)) { # 75 percent sig
				segments(time[i], d.ll75[i], time[i], d.ul75[i], lwd = 5, col = "grey30")
			}
			# Actual response
			points(time, d.central, lwd = 4)	
		}
	} else if(response == "shock.effect.decay") { # The decay in the shock is the absolute value, period by period, of diffs (which die off but could oscillate)
		for(i in 1:length(z$d.central)) { # correct for positive/negative
			if(z$d.central[i] < 0) {# if it's negative
				z$d.central[i] <- z$d.central[i] * (-1)
				temp <- c(z$d.ll95[i], z$d.ll90[i], z$d.ll75[i], z$d.ul75[i], z$d.ul90[i], z$d.ul95[i])
				temp <- temp * (-1)
				z$d.ul95[i] <- temp[1] # And reverse the order
				z$d.ul90[i] <- temp[2]
				z$d.ul75[i] <- temp[3]
				z$d.ll75[i] <- temp[4]
				z$d.ll90[i] <- temp[5]
				z$d.ll95[i] <- temp[6]
			}
		}
		plot(z$time, z$d.ll95, type = "n",
			ylim = ifelse(c(is.null(ylim), is.null(ylim)), c(min(z$d.ll95), max(z$d.ul95)), ylim),
			ylab = ifelse(is.null(ylab), "Shock Effect Value", ylab), xlab = ifelse(is.null(xlab), "Time", xlab), ...)
		if(type == "area") { 
			polygon(c(z$time, rev(z$time)), c(z$d.ul95, rev(z$d.ll95)), col = ifelse(bw == FALSE, "skyblue1", "grey70"), border = NA) # 95
			polygon(c(z$time, rev(z$time)), c(z$d.ul90, rev(z$d.ll90)), col = ifelse(bw == FALSE, "skyblue3", "grey50"), border = NA) # 90
			polygon(c(z$time, rev(z$time)), c(z$d.ul75, rev(z$d.ll75)), col = "grey30", border = NA) # 75
			# Actual response
			lines(z$time, z$d.central, lty = 2, lwd = 3)
		} else { # if it's a spikeplot
			for(i in 1:length(z$time)) { # 95 percent sig
				segments(z$time[i], z$d.ll95[i], z$time[i], z$d.ul95[i], lwd = 1, col = ifelse(bw == FALSE, "skyblue1", "grey70"))
			}
			for(i in 1:length(z$time)) { # 90 percent sig
				segments(z$time[i], z$d.ll90[i], z$time[i], z$d.ul90[i], lwd = 3, col = ifelse(bw == FALSE, "skyblue3", "grey50"))
			}
			for(i in 1:length(z$time)) { # 75 percent sig
				segments(z$time[i], z$d.ll75[i], z$time[i], z$d.ul75[i], lwd = 5, col = "grey30")
			}
			# Actual response
			points(z$time, z$d.central, lwd = 4)	
		}
	}
	else if(response == "cumulative.abs.diffs") { # If it's cumulative differences in Y
		if(identical(x$rawsims, NULL)) {
			stop("dynardl object must have fullsims = TRUE to track cumulative absolute differences.")
		}
		a <- x$rawsims[,start.period:ncol(x$rawsims)]
		# Here, diffs sims is going to be 2:length of simulation, with the first period NA. This is to sync it with the d.central output
		diff.sims <- matrix(rep(NA, nrow(a)*(ncol(a) - 1)), nrow = nrow(a)) # Last column is central tendency
		# Track the differences for each simulation
		for(i in 2:ncol(diff.sims)) {
			diff.sims[,i] <- a[,i] - a[,(i - 1)]
		}
		is.changing.test <- is.changing
		if(identical(is.changing, NULL)) {
			is.changing.test <- 1
		}		
		temp.ll95 <- temp.ll90 <- temp.ll75 <- temp.ul75 <- temp.ul90 <- temp.ul95 <- temp.central <- rep(NA, ncol(diff.sims))
		d.ll95 <- d.ul95 <- d.ll90 <- d.ul90 <- d.ll75 <- d.ul75 <- d.central <- rep(NA, ncol(diff.sims))		
		for(i in 2:ncol(diff.sims)) {
			if((i + start.period - 1) < x$simulation$shocktime[1]) { # If it's in the equilibriating period before the shock
				temp.ll95[i] <- quantile(diff.sims[,i], 0.025, na.rm = T) # Preserve the regular diffs, NOT absolute, since, they're noise
				temp.ll90[i] <- quantile(diff.sims[,i], 0.050, na.rm = T)
				temp.ll75[i] <- quantile(diff.sims[,i], 0.125, na.rm = T)
				temp.ul75[i] <- quantile(diff.sims[,i], 0.875, na.rm = T)
				temp.ul90[i] <- quantile(diff.sims[,i], 0.950, na.rm = T)
				temp.ul95[i] <- quantile(diff.sims[,i], 0.975, na.rm = T)
				temp.central[i] <- ifelse(identical(x$rawsims$central[1], "median"), median(diff.sims[,i], na.rm = T), mean(diff.sims[,i], na.rm = T))
				d.central[i] <- sum(temp.central[1:i], na.rm = T)
			} else {
				if((is.changing.test - start.period + 1) >= i) { # if it is moving
					if(mean(diff.sims[,i]) > 0) { # if the movement is positive
						temp.central[i] <- ifelse(identical(x$rawsims$central[1], "median"), median(diff.sims[,i], na.rm = T), mean(diff.sims[,i], na.rm = T))
						d.central[i] <- sum(temp.central[1:i], na.rm = T)
						if(abs.errors == "cumulate") {
							temp.ll95[i] <- temp.ll95[i-1] + (quantile(diff.sims[,i], 0.025, na.rm = T))
							temp.ll90[i] <- temp.ll90[i-1] + (quantile(diff.sims[,i], 0.050, na.rm = T))
							temp.ll75[i] <- temp.ll75[i-1] + (quantile(diff.sims[,i], 0.125, na.rm = T))
							temp.ul75[i] <- temp.ul75[i-1] + (quantile(diff.sims[,i], 0.875, na.rm = T))
							temp.ul90[i] <- temp.ul90[i-1] + (quantile(diff.sims[,i], 0.950, na.rm = T))
							temp.ul95[i] <- temp.ul95[i-1] + (quantile(diff.sims[,i], 0.975, na.rm = T))
						} else { # if it's none or within, do the same for the plot scale
							temp.ll95[i] <- d.central[i] - abs(temp.central[i] - quantile(diff.sims[,i], 0.025, na.rm = T))
							temp.ll90[i] <- d.central[i] - abs(temp.central[i] - quantile(diff.sims[,i], 0.050, na.rm = T))
							temp.ll75[i] <- d.central[i] - abs(temp.central[i] - quantile(diff.sims[,i], 0.125, na.rm = T))
							temp.ul75[i] <- d.central[i] + (quantile(diff.sims[,i], 0.875, na.rm = T) - temp.central[i])
							temp.ul90[i] <- d.central[i] + (quantile(diff.sims[,i], 0.950, na.rm = T) - temp.central[i])
							temp.ul95[i] <- d.central[i] + (quantile(diff.sims[,i], 0.975, na.rm = T) - temp.central[i])
						}
					} else { # if the movement is negative
						temp.central[i] <- ifelse(identical(x$rawsims$central[1], "median"), median(diff.sims[,i], na.rm = T)*(-1), mean(diff.sims[,i], na.rm = T)*(-1))
						d.central[i] <- sum(temp.central[1:i], na.rm = T)
						if(abs.errors == "cumulate") {
							temp.ul95[i] <- temp.ul95[i-1] + (quantile(diff.sims[,i], 0.025, na.rm = T)*(-1))
							temp.ul90[i] <- temp.ul90[i-1] + (quantile(diff.sims[,i], 0.050, na.rm = T)*(-1))
							temp.ul75[i] <- temp.ul75[i-1] + (quantile(diff.sims[,i], 0.125, na.rm = T)*(-1))
							temp.ll75[i] <- temp.ll75[i-1] + (quantile(diff.sims[,i], 0.875, na.rm = T)*(-1))
							temp.ll90[i] <- temp.ll90[i-1] + (quantile(diff.sims[,i], 0.950, na.rm = T)*(-1))
							temp.ll95[i] <- temp.ll95[i-1] + (quantile(diff.sims[,i], 0.975, na.rm = T)*(-1))
						} else { # if it's none or within, do the same for the plot scale						
							temp.ul95[i] <- d.central[i] + (quantile(diff.sims[,i], 0.025, na.rm = T)*(-1) - temp.central[i])
							temp.ul90[i] <- d.central[i] + (quantile(diff.sims[,i], 0.050, na.rm = T)*(-1) - temp.central[i])
							temp.ul75[i] <- d.central[i] + (quantile(diff.sims[,i], 0.125, na.rm = T)*(-1) - temp.central[i])
							temp.ll75[i] <- d.central[i] - abs(temp.central[i] - quantile(diff.sims[,i], 0.875, na.rm = T)*(-1))
							temp.ll90[i] <- d.central[i] - abs(temp.central[i] - quantile(diff.sims[,i], 0.950, na.rm = T)*(-1))
							temp.ll95[i] <- d.central[i] - abs(temp.central[i] - quantile(diff.sims[,i], 0.975, na.rm = T)*(-1))
						}
					}
				} else { # if it's not moving, by tol or last.period
						temp.ll95[i] <- temp.ll95[i-1]
						temp.ll90[i] <- temp.ll90[i-1]
						temp.ll75[i] <- temp.ll75[i-1]
						temp.ul75[i] <- temp.ul75[i-1]
						temp.ul90[i] <- temp.ul90[i-1]
						temp.ul95[i] <- temp.ul95[i-1]
						temp.central[i] <- 0
						d.central[i] <- sum(temp.central[1:i], na.rm = T)
					}
				}
			d.ll95[i] <- temp.ll95[i] 
			d.ll90[i] <- temp.ll90[i] 
			d.ll75[i] <- temp.ll75[i] 
			d.ul75[i] <- temp.ul75[i] 
			d.ul90[i] <- temp.ul90[i] 
			d.ul95[i] <- temp.ul95[i] 
		}
		time <- seq(start.period, (ncol(diff.sims) + start.period - 1), 1)
		plot(time, d.ll95, type = "n", 
			ylim = ifelse(c(is.null(ylim), is.null(ylim)), c(min(d.ll95, na.rm = T), max(d.ul95, na.rm = T)), ylim),
			ylab = ifelse(is.null(ylab), "Cumulative Absolute Change in Y Value", ylab), xlab = ifelse(is.null(xlab), "Time", xlab))
		if(type == "area") {
			if(abs.errors != "none") {
				polygon(c(time, rev(time)), c(d.ul95, rev(d.ll95)), col = ifelse(bw == FALSE, "skyblue1", "grey70"), border = NA) # 95
				polygon(c(time, rev(time)), c(d.ul90, rev(d.ll90)), col = ifelse(bw == FALSE, "skyblue3", "grey50"), border = NA) # 90
				polygon(c(time, rev(time)), c(d.ul75, rev(d.ll75)), col = "grey30", border = NA) # 75
			}
			# Actual response
			lines(time, d.central, lty = 2, lwd = 3)
		} else { # if it's a spikeplot
			if(abs.errors != "none") {
				for(i in 1:length(time)) { # 95 percent sig
					segments(time[i], d.ll95[i], time[i], d.ul95[i], lwd = 1, col = ifelse(bw == FALSE, "skyblue1", "grey70"))
				}
				for(i in 1:length(time)) { # 90 percent sig
					segments(time[i], d.ll90[i], time[i], d.ul90[i], lwd = 3, col = ifelse(bw == FALSE, "skyblue3", "grey50"))
				}
				for(i in 1:length(time)) { # 75 percent sig
					segments(time[i], d.ll75[i], time[i], d.ul75[i], lwd = 5, col = "grey30")
				}
			}
			# Actual response
			points(time, d.central, lwd = 4)
		}
		if(!(identical(last.period, NULL))) {
			warning(paste(paste("Cumulative absolute effects assumed to be noise at t = "), paste(is.changing), paste(" (last.period)."), sep = ""))		
		} else {
			warning(paste(paste("Cumulative absolute effects assumed to be noise (by tolerance) at t = "), paste(is.changing.test), paste("."), sep = ""))
		}
		if(identical(is.changing, NULL)) { # If Y never responds
			warning("Y does not move beyond the tolerance in the simulation. Reconsider the tolerance, or investigate if Y responds to the shockvar in the dynardl model.")
		} else if(is.changing == length(z$d.central)) { # If the simulation isn't long enough, potentially
			warning("Y might still be changing (has not met tolerance) at the end of the simulation. Consider lengthening the simulation in dynardl or adjusting the tolerance.")
		}
	}
}


#############################################
# ---------(11) dynardl.all.plots ----------#
#############################################
#' Combine all of the potential plots of a simulated response in a \code{\link{dynardl}} model
#' @param x a \code{dynardl} model with a simulation to be plotted. Since all plots includes absolute cumulative differences, \code{fullsims} must be \code{TRUE} in the \code{dynardl} simulation
#' @param type whether the plot should be an area plot (\code{area}) or a spike plot (\code{spike})
#' @param bw should the colors be in black and white (for publication)? The default is \code{FALSE}
#' @param tol when deciding when to stop calculating the absolute value of the shocks to the dependent variable, you can specify the minimum amount of movement required to qualify as a non-noise change over time periods (for calculating absolute cumulative differences). The default is 0.1 percent of the mean of the dependent variable. Specify a \code{tol} or a \code{last.period}. If both are specified, \code{last.period} overrides \code{tol}
#' @param start.period which period of the simulation to begin the plot with. You can view the equilibriating behavior of the dependent variable, or you can skip forward in time (maybe to just before the shock). The default is \code{1} (the first period of the simulation)
#' @param last.period when deciding when to stop calculating the absolute value of the shocks to the dependent variable, you can specify a specific period in which to stop calculating absolute cumulative differences. Specify a \code{tol} or a \code{last.period}. If both are specified, \code{last.period} overrides \code{tol}
#' @param abs.errors when calculating confidence for the absolute cumulative effect, should differences accumulate in each time time period (\code{cumulate}, which could be explosive if the error in the model is large), should differences be observed at each time (\code{within.period}, which will have smaller values in equilibrium than when changing), or should only the values be plotted (\code{none})
#' @param ylab a user-defined y-label to be used instead of the default (use caution, as it will be passed to all plots)
#' @param xlab a user-defined x-label to be used instead of the default (use caution, as it will be passed to all plots)
#' @param ylim a user-defined y-limit to be used instead of the default (for instance, for shared axes. Use caution, as it will be passed to all plots)
#' @param ... other arguments to be passed to the call to plot. Use caution, as they will be passed to all plots
#' @return a 2 x 3 grid of the plots of the simulated dynardl model effects plots
#' @details
#' When running \code{dynardl}, \code{simulate} must be \code{TRUE} so that there is a simulation to plot. Also, \code{fullsims} must be \code{TRUE} as the plot will contain absolute cumulative differences. See \code{\link{dynardl.simulation.plot}} for arguments to the individual plotting types
#' @importFrom graphics lines plot points polygon segments par
#' @author Soren Jordan and Andrew Q. Philips
#' @keywords utilities
#' @examples
#' # Using the ineq data in dynamac
#' # Shocking Income Top 10
#' set.seed(1)
#' ardl.model <- dynardl(concern ~ incshare10 + urate, data = ineq, 
#'        lags = list("concern" = 1, "incshare10" = 1),
#'        diffs = c("incshare10", "urate"), 
#'        lagdiffs = list("concern" = 1),
#'        ec = TRUE, simulate = TRUE, range = 30,
#'        shockvar = "incshare10", fullsims = TRUE)
#' 
#' # Shows all of the potential responses
#' dynardl.all.plots(ardl.model)	
#' # Same plot, but with spikeplot
#' dynardl.all.plots(ardl.model, type = "spike")  
#' # Grayscale plots
#' dynardl.all.plots(ardl.model, bw = TRUE)	 
#' @export

dynardl.all.plots <- function(x, type = "area", bw = FALSE, last.period = NULL, start.period = 1, tol = (abs(x$model$ymean) * 0.01), abs.errors = "none", ylim = NULL, xlab = NULL, ylab = NULL, ...) {
	if(!(identical(class(x), "dynardl"))) { # Check if data are from dynardl
		stop("To plot simulation, provide a dynardl model object (with a simulation).")
	}
	if(x$model$simulate == FALSE) {
		stop("dynardl object does not include simulation to plot.")
	}
	if(identical(x$rawsims, NULL)) {
		stop("dynardl object must have fullsims = TRUE to track cumulative absolute differences.")
	}
	par(mfrow = c(2, 3))
	dynardl.simulation.plot(x, response = "levels", type = type, bw = bw, tol = tol, last.period = last.period, start.period = start.period, ylim = ylim, ylab = ylab, xlab = xlab, ...)
	dynardl.simulation.plot(x, response = "levels.from.mean", type = type, bw = bw, tol = tol, last.period = last.period, start.period = start.period, ylim = ylim, ylab = ylab, xlab = xlab,  ...)
	dynardl.simulation.plot(x, response = "diffs", type = type, bw = bw, tol = tol, last.period = last.period, start.period = start.period, ylim = ylim, ylab = ylab, xlab = xlab, ...)
	dynardl.simulation.plot(x, response = "shock.effect.decay", type = type, bw = bw, tol = tol, last.period = last.period, start.period = start.period, ylim = ylim, ylab = ylab, xlab = xlab, ...)
	dynardl.simulation.plot(x, response = "cumulative.diffs", type = type, bw = bw, tol = tol, last.period = last.period, start.period = start.period, ylim = ylim, ylab = ylab, xlab = xlab, ...)
	dynardl.simulation.plot(x, response = "cumulative.abs.diffs", type = type, bw = bw, tol = tol, last.period = last.period, start.period = start.period, abs.errors = abs.errors, ylim = ylim, ylab = ylab, xlab = xlab, ...)
}


###################################
# ---(12) dynardl.totaleffect ----#		# HAS NOT BEEN FIXED FOR NEW CHANGE TO CUM. ABS.DIFFS; no warning for LDV, etc. etc. 
###################################
dynardl.totaleffect <- function(x, last.period = NULL, tol = (abs(x$model$ymean) * 0.001), round.to = 3, object.out = FALSE) {
	if(x$model$simulate == FALSE) {
		stop("dynardl object does not include simulation to calculate.")
	}
	if(identical(x$rawsims, NULL)) {
		stop("dynardl object must have fullsims = TRUE to calculate total effects.")
	}
	if(!(identical(last.period, NULL))) { # If we're picking a time period to stop
		if(last.period > length(x$simulation$central)) {
			warning(paste(paste("last.period requested exceeds simulation range. Calculating on"), paste(length(x$simulation$central)), paste("periods."), sep = " "))
			is.changing <- last.period <- length(x$simulation$central)
		}
		warning(paste(paste("Cumulative absolute change in Y calculated on "), paste(last.period), paste(" periods: regardless of if Y might still be moving OR calculation might include noise in cumulative effect."), sep = ""))
		is.changing <- last.period
	} else { # If we're calculating based on tolerance
		is.changing <- NULL
		for(i in 1:length(x$simulation$d.central)) { # Test if it's changing
			if(abs(x$simulation$d.central[i]) > tol) {
				is.changing <- i
			}
		}
	}
	# First, the changes from the mean (to see the final movement offset by the absolutes)
	from.mean.ll95 <- x$simulation$ll95 - x$model$ymean
	from.mean.ul95 <- x$simulation$ul95 - x$model$ymean
	from.mean.ll90 <- x$simulation$ll90 - x$model$ymean
	from.mean.ul90 <- x$simulation$ul90 - x$model$ymean
	from.mean.ll75 <- x$simulation$ll75 - x$model$ymean
	from.mean.ul75 <- x$simulation$ul75 - x$model$ymean
	from.mean.central <- x$simulation$central - x$model$ymean
	# Now with the same calculations as abs code in graphics
	diff.sims <- cum.diff.sims <- temp.abs.diff.sims <- cum.abs.diff.sims <- matrix(rep(NA, nrow(x$rawsims)*(ncol(x$rawsims) - 1)), nrow = nrow(x$rawsims)) # Last column is central tendency
	for(i in 2:ncol(diff.sims)) {
		diff.sims[,i] <- x$rawsims[,i] - x$rawsims[,(i - 1)]
	}
	# First period: no cumulative difference
	# Up until the shocktime, we're going to take regular diffs, NOT absolute, since, they're noise	
	temp.abs.diff.sims[,2] <- cum.abs.diff.sims[,2] <- diff.sims[,2]
	for(i in 3:ncol(temp.abs.diff.sims)) {
		if(i < x$simulation$shocktime[1]) { # If it's in the equilibriating period before the shock
			temp.abs.diff.sims[,i] <- diff.sims[,i] # Preserve the regular diffs, NOT absolute, since, they're noise
		} else {
			if(is.changing >= i) { # If it is still moving on the differences from the beginning,
				temp.abs.diff.sims[,i] <- abs(diff.sims[,i]) # Preserve the ABSOLUTE changes for that time period
			} else { # If the changes aren't `real', meaning below our tolerance
				temp.abs.diff.sims[,i] <- 0 # Replace the `change' (noise) with nothing
			}
		}
		# Now: the sims we're going to keep and graph: sum over all of them to now
		cum.abs.diff.sims[,i] <- rowSums(temp.abs.diff.sims[,1:i], na.rm = T)
	}
	# Output
	cum.abs.ll95 <- quantile(cum.abs.diff.sims[,is.changing], 0.025, na.rm = T)
	cum.abs.ll90 <- quantile(cum.abs.diff.sims[,is.changing], 0.050, na.rm = T)
	cum.abs.ll75 <- quantile(cum.abs.diff.sims[,is.changing], 0.125, na.rm = T)
	cum.abs.ul75 <- quantile(cum.abs.diff.sims[,is.changing], 0.875, na.rm = T)
	cum.abs.ul90 <- quantile(cum.abs.diff.sims[,is.changing], 0.950, na.rm = T)
	cum.abs.ul95 <- quantile(cum.abs.diff.sims[,is.changing], 0.975, na.rm = T)
	cum.abs.central <- ifelse(identical(x$rawsims$central[1], "median"), median(cum.abs.diff.sims[,is.changing], na.rm = T), mean(cum.abs.diff.sims[,is.changing], na.rm = T))
	flush.console()
	if(!(identical(last.period, NULL))) {
		change.message <- paste(paste("Total effects calculated on "), paste(last.period), paste(" periods."), sep = "")
	} else {
		change.message <- ifelse(is.changing == ncol(temp.abs.diff.sims), # If the simulation isn't long enough, potentially
				paste("Differences in Y might still be changing (have not met tolerance) at the end of the simulation. Consider lengthening the simulation in dynardl or adjusting the tolerance."),
				paste(paste("Differences in Y met the tolerance (stopped moving) between the "), paste(is.changing - 1), paste(" and "), paste(is.changing), paste(" time periods (including pre-shock)."), sep = ""))
	}
	cat("\n",
		"------------------------------------------------------", "\n",
		paste(change.message), "\n",
		"------------------------------------------------------", "\n",
		paste(paste("Ending level of Y (away from mean): "), paste(round(from.mean.central[is.changing], digits = round.to)), sep = ""), "\n",
		paste(paste("95% percentile of ending level of Y (away from mean): ["), paste(round(from.mean.ll95[is.changing], digits = round.to)), paste(", "), paste(round(from.mean.ul95[is.changing], digits = round.to)), paste("]."), sep = ""), "\n",
		paste(paste("90% percentile of ending level of Y (away from mean): ["), paste(round(from.mean.ll90[is.changing], digits = round.to)), paste(", "), paste(round(from.mean.ul90[is.changing], digits = round.to)), paste("]."), sep = ""), "\n",
		paste(paste("75% percentile of ending level of Y (away from mean): ["), paste(round(from.mean.ll75[is.changing], digits = round.to)), paste(", "), paste(round(from.mean.ul75[is.changing], digits = round.to)), paste("]."), sep = ""), "\n",
		"------------------------------------------------------", "\n",
		paste(paste("Ending cumulative absolute change (total movement in Y): "), paste(round(cum.abs.central, digits = round.to)), paste("."), sep = ""), "\n",
		paste(paste("95% percentile of ending cumulative absolute changes: ["), paste(round(cum.abs.ll95, digits = round.to)), paste(", "), paste(round(cum.abs.ul95, digits = round.to)), paste("]."), sep = ""), "\n",
		paste(paste("90% percentile of ending cumulative absolute changes: ["), paste(round(cum.abs.ll90, digits = round.to)), paste(", "), paste(round(cum.abs.ul90, digits = round.to)), paste("]."), sep = ""), "\n",
		paste(paste("75% percentile of ending cumulative absolute changes: ["), paste(round(cum.abs.ll75, digits = round.to)), paste(", "), paste(round(cum.abs.ul75, digits = round.to)), paste("]."), sep = ""), "\n",
		"------------------------------------------------------", "\n",
		paste(paste("Calculated on "), paste(is.changing - x$simulation$shocktime[1]), paste(" periods after the shock (at t = "), paste(x$simulation$shocktime[1]), paste(")."), sep = ""), "\n",
		"------------------------------------------------------", "\n")
	# Now issue a few warnings
	if(identical(is.changing, NULL)) { # If Y never responds
		warning("Differences in Y do not move beyond the tolerance in the simulation. Reconsider the tolerance, or investigate if Y responds to the shockvar in the dynardl model.")
	} else if(is.changing == ncol(temp.abs.diff.sims) & !(identical(last.period, NULL))) { # If the simulation isn't long enough, potentially
		warning("Differences in Y might still be changing (have not met tolerance) at the end of the simulation. Consider lengthening the simulation in dynardl or adjusting the tolerance.")
	}
	if(object.out == TRUE) {
		level.from.mean <- matrix(c(from.mean.central[is.changing], 
						from.mean.ll95[is.changing], from.mean.ul95[is.changing], 
						from.mean.ll90[is.changing], from.mean.ul90[is.changing],
						from.mean.ll75[is.changing], from.mean.ul75[is.changing]), nrow = 1)
		cumulative.abs.changes <- matrix(c(cum.abs.central, 
						cum.abs.ll95, cum.abs.ul95, 
						cum.abs.ll90, cum.abs.ul90,
						cum.abs.ll75, cum.abs.ul75), nrow = 1)
		colnames(level.from.mean) <- colnames(cumulative.abs.changes) <- c("central", "ll95", "ul95", "ll90", "ul90", "ll75", "ul75")
		out <- structure(list(level.quantities = level.from.mean, cumulative.quantities = cumulative.abs.changes))
		out
	}	
}	






















































































































###################################
# -----(XX) dynardl.effects ------#
###################################

dynardl.effects <- function(x, tol = 0.025, period = NULL, x.lag = NULL, object.out = FALSE) {
	if(x$model$simulate == FALSE) {
		stop("dynardl object does not include simulation to plot.")
	}
	shocktime <- x$simulation$shocktime[1]
	changes <- abs(dshift(x$simulation$mean))
	dists <- abs.dists <- rep(NA, length(changes))
	innov <- abs.innov <- rep(NA, length(changes))
	triangles <- abs.triangles <- rep(NA, length(changes))
	triangles.rectangles <- abs.triangles.rectangles <- rep(NA, length(changes))
	sig.changes <- rep(NA, length(changes))
	demeaned <- x$simulation$mean - x$model$ymean
	for(i in 2:length(sig.changes)) {
		sig.changes[i] <- ifelse(changes[i] >= tol*(max(x$simulation$mean)-min(x$simulation$mean)), TRUE, FALSE) # is the new change more or less than (tolerance)% of the total?
	}
	for(i in 2:length(triangles)) {
		if((demeaned[i] > 0) & (demeaned[i-1] > 0)) { # If both Y values are positive
			abs.dists[i] <- dists[i] <- demeaned[i]	 # Vertical measure: just the value of the demeaned series (the spike)
			abs.innov[i] <- innov[i] <- changes[i] # Unit over unit change in the spike level
			abs.triangles[i] <- triangles[i] <- changes[i]*0.5 # 1/2*base(always 1 unit) * height: change in value t-1 to t
			abs.triangles.rectangles[i] <- triangles.rectangles[i] <- triangles[i] + min(demeaned[i], demeaned[i-1]) # It's the trianlge plus the rectangle underneath it: l x w (w = 1, height is the lesser of the two verticals from the mean)
		} else {
			if((demeaned[i] < 0) & (demeaned[i-1] < 0)) {	# If both Y values are negative
				dists[i] <- demeaned[i] # Vertical measure: just the value of the demeaned series (the spike)
				abs.dists[i] <- abs(dists[i]) # Don't make it negative, as it's counting ``total effect''
				innov[i] <- (-1)*changes[i] # Negative change from last period
				abs.innov[i] <- abs(innov[i]) # Total change from previous 
				triangles[i] <- (-1)*changes[i]*0.5 # 1/2*base(always 1 unit) * height
				abs.triangles[i] <- abs(triangles[i])
				triangles.rectangles[i] <- triangles[i] + max(demeaned[i], demeaned[i-1])	 # The triangle plus the rectangle underneath which is the GREATER (less negative) of the two	
				abs.triangles.rectangles[i] <- abs(triangles.rectangles[i])		
			} else {
				# X runs one unit, total change is changes[i], so slope of full line is changes[i]
				#  Chunk up the line by reducing the slope into the parts of changes that are above and below changes
				if(demeaned[i] < 0) { # If the second observation is below zero
					cross.over <- (demeaned[i-1]/changes[i])
					pos.triangle <- demeaned[i-1]*cross.over*0.5
					neg.triangle <- demeaned[i]*(1-cross.over)*0.5
					innov[i] <- dists[i] <- demeaned[i-1] + demeaned[i] # The net (positive/negative) of the two spikes
					abs.innov[i] <- abs.dists[i] <- demeaned[i-1] + abs(demeaned[i]) # Make the second observation positive!
					triangles[i] <- triangles.rectangles[i] <- pos.triangle + neg.triangle # No rectangles below: crossover
					abs.triangles[i] <- abs.triangles.rectangles[i] <- pos.triangle + abs(neg.triangle)
				} else {
					cross.over <- (demeaned[i]/changes[i])
					pos.triangle <- demeaned[i]*cross.over*0.5
					neg.triangle <- demeaned[i-1]*(1-cross.over)*0.5
					innov[i] <- dists[i] <- demeaned[i-1] + demeaned[i] # The net (positive/negative) of the two spikes
					abs.innov[i] <- abs.dists[i] <- abs(demeaned[i-1]) + demeaned[i] # Make the first observation positive!
					triangles[i] <- triangles.rectangles[i] <- pos.triangle + neg.triangle # No rectangles below: crossover
					abs.triangles[i] <- abs.triangles.rectangles[i] <- pos.triangle + abs(neg.triangle)	
				}
			}
		}
	}
	areas <- data.frame(x$simulation$mean, demeaned, dists, abs.dists, innov, abs.innov, triangles, abs.triangles, triangles.rectangles, abs.triangles.rectangles)
	names(areas) <- c("Y.sim.mean", "Y.sim.demeaned", "spike.lengths", "abs.spike.lengths", "innovations", "abs.innovations",
		"innovation.area", "abs.innovation.area", "total.area", "abs.total.area")
	t.spikes <- a.t.spikes <- t.innov <- a.t.innov <- t.trian <- a.t.trian <- t.trirect <- a.t.trirect <- NULL
	# Now actually calculate 
	# When's the last significant change (according to tolerance?). Used for testing
	which.to.add <- seq(1, length(demeaned), 1)
	last.sign.period <- max(subset(which.to.add, sig.changes == TRUE))

	if(is.null(period)) { # If we're not calculating based on a period
		last.period <- last.sign.period
		t.spikes <- sum(dists[shocktime:last.period])
		a.t.spikes <- sum(abs.dists[shocktime:last.period])
		t.innov <- sum(innov[shocktime:last.period])
		a.t.innov <- sum(abs.innov[shocktime:last.period])		
		t.trian <- sum(triangles[shocktime:last.period])
		a.t.trian <- sum(abs.triangles[shocktime:last.period])
		t.trirect <- sum(triangles.rectangles[shocktime:last.period])
		a.t.trirect <- sum(abs.triangles.rectangles[shocktime:last.period])
	} else { # If we ARE calculating based on a period
		last.period <- shocktime + period
		if(last.period > max(x$simulation$time)) { # If the user wants more periods than are available ...
			last.period <- max(x$simulation$time)
			warning(paste(paste("Length of effect (period) requested exceeds simulation length. Effect calculated for"), paste(last.period - shocktime), paste("periods."), sep = " "))	
		}
		if(last.period < last.sign.period) { # If the user wants a period but Y could still be responding ...
			warning(paste(paste("Y might still be changing after"), paste(period), paste("period(s). Consider lengthening the number of periods included in calculated effect."), sep = " "))	
		}
		t.spikes <- sum(dists[shocktime:last.period])
		a.t.spikes <- sum(abs.dists[shocktime:last.period])
		t.innov <- sum(innov[shocktime:last.period])
		a.t.innov <- sum(abs.innov[shocktime:last.period])
		t.trian <- sum(triangles[shocktime:last.period])
		a.t.trian <- sum(abs.triangles[shocktime:last.period])
		t.trirect <- sum(triangles.rectangles[shocktime:last.period])
		a.t.trirect <- sum(abs.triangles.rectangles[shocktime:last.period])
	}
	# Now issue a few warnings
	if(last.sign.period == max(which.to.add)) { # If the simulation isn't long enough, potentially
		warning(paste(paste("Y might still be changing at the end of the simulation. Consider lengthening the simulation beyond range ="), paste(max(which.to.add)), paste("periods in dynardl."), sep = " "))
	}
	if(!(min(x$simulation$ul75) < max(x$simulation$ll75))) { # If there is overlap in the significant regions ...
		warning(paste("Use caution, as the effect of X on Y might not be meaningfully different from 0: upper 75% interval overlaps with lower 75% interval."))
	}
	if(is.null(x.lag)) {
		immediate <- demeaned[shocktime] - demeaned[shocktime - 1]
	} else {
		immediate <- demeaned[shocktime + x.lag] - demeaned[shocktime + x.lag - 1]
	}
	quantities <- data.frame(immediate, t.spikes, a.t.spikes, t.innov, a.t.innov, t.trian, a.t.trian, t.trirect, a.t.trirect)
	names(quantities) <- c("immediate.effect", "sum.spike.lengths", "sum.abs.spike.lengths", "sum.innovations", "sum.abs.innovations",
		"sum.innovation.area", "sum.abs.innovation.area", "sum.total.area", "sum.abs.total.area")
	out <- list(areas, quantities)
	names(out) <- c("areas", "quantities")
	
	# Output
	flush.console()
	cat("\n",
		"------------------------------------------------------", "\n")
		if(is.null(x.lag)) {
			cat(paste(" Immediate (shock period) effect:", round(out$quantities$immediate.effect, digits = 3), sep = " "), "\n")
		} else {
			cat(paste(" Immediate (shock period) effect:", round(out$quantities$immediate.effect, digits = 3), "adjusted by", x.lag, "period(s) for the lag of X.", sep = " "), "\n")
		}
		cat(
		" ------------------------------------------------------", "\n",
		paste("Net sum of spikes:", round(out$quantities$sum.spike.lengths, digits = 3), sep = " "), "\n",
		paste("Absolute sum of spikes:", round(out$quantities$sum.abs.spike.lengths, digits = 3), sep = " "), "\n",
		paste("Net sum of innovations:", round(out$quantities$sum.innovations, digits = 3), sep = " "), "\n",
		paste("Absolute sum of innovations:", round(out$quantities$sum.abs.innovations, digits = 3), sep = " "), "\n",
		paste("Net area of innovations:", round(out$quantities$sum.innovation.area, digits = 3), sep = " "), "\n",
		paste("Absolute area of innovations:", round(out$quantities$sum.abs.innovation.area, digits = 3), sep = " "), "\n",
		paste("Net sum of area under curve:", round(out$quantities$sum.total.area, digits = 3), sep = " "), "\n",
		paste("Absolute sum of area under curve:", round(out$quantities$sum.abs.total.area, digits = 3), sep = " "), "\n",
		paste("Calculated on", paste(last.period - shocktime), paste("periods."), sep = " "), "\n",
		"------------------------------------------------------", "\n")
	if(object.out == TRUE) {
		out
	}	
}