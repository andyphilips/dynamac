--- 
title: 'Exploring meaningful visual effects and quantities of interest from dynamic models through dynamac'  
authors:
  - name: Soren Jordan 
    affiliation: 1 
    orcid: 0000-0003-4201-1085 
  - name: Andrew Q. Philips 
    affiliation: 2 
affiliations: 
  - index: 1 
    name: Assistant Professor, Department of Political Science, Auburn University 
  - index: 2 
    name: Assistant Professor, Department of Political Science, University of Colorado Boulder 
date: 8 October 2020 
tags: 
  - time series
  - methodology
  - long-run effects
  - simulations
  - quantities of interest
bibliography: jordan-philips.bib 
---

# Summary 

dynamac [-@dynamacCRAN] implements a framework of estimating and interpreting autoregressive distributed lag (ARDL) models in R. dynamac uses stochastic simulation techniques [@jordan2018cointegration;@jordan2018dynamic] to easily recover traditional quantities of interest, such as short- and long-run effects, even from complex dynamic specifications, including models with cointegration. These simulation techniques bring a wider set of inferences from complex models to users in fields as diverse as environmental science [@jcp] and economics [@te]. Users can make dynamic inferences from their models, including measures of uncertainty, without the need for any formulae or algebraic solutions. Although other packages may estimate ARDL models (e.g., dynlm [@dynlmCRAN], which only provides model estimates; ardl [@ardlCRAN], which uses dynlm for estimation but provides support for the bounds test for cointegration and automated lag selection; dlagM [@dlagMCRAN], which allows for forecasting), dynamac is unique in that it provides post-estimation diagnostics for autocorrelation in the residuals and is designed specifically for providing inferences through counterfactual simulation.

The procedure first estimates the parameters using an ARDL model:
\begin{align}\label{eq:ardlgeneral}
\begin{split}
	y_t,\Delta y_t =  &\alpha_0 + \delta T + \sum_{p=1}^P \phi_p y_{t-p} + \sum_{l_1 = 0}^{L_1} \theta_{1l}x_{1,t-l_1} + \cdots + \sum_{l_k = 0}^{L_k} \theta_{kl}x_{t-l_k} + \\
	& \sum_{m=1}^M \alpha_m \Delta y_{t-m} + \sum_{q_1=0}^{Q_1} \beta_{1q_{1}} \Delta x_{1,t-q{_1}} + \cdots + \sum_{q_k=0}^{Q_k} \beta_{kq_k} \Delta x_{k,t-q{_k}} + \epsilon_t
\end{split}
\end{align} The dependent variable (appearing in level form or first-differences) is a function of a constant, a deterministic linear trend $T$, up to $P$ and $L$ lags of the dependent and independent variables, up to $M$ and $Q$ lags of the first-differenced dependent and independent variables, and  error term $\epsilon_t$.  Users should use a combination of theory, information criteria, unit-root, cointegration and residual-based tests to arrive at a more restrictive specification of Equation \ref{eq:ardlgeneral}.

Next, through `dynardl()`, the estimated parameters are used to draw $s$ simulations from a multivariate normal distribution, with mean $\hat{\beta}$ and variance obtained from the estimated variance-covariance matrix. Stable predicted values for $y$ using the $s$ simulations of $\hat{\beta}$ are then created, using starting values of the independent variables (usually means for continuous variables or modes for categorical ones). Users can choose from either expected values or---by incorporating fundamental model uncertainty---predicted values. Next, at a defined time period $t$, one of the independent variables is "shocked" by an amount determined by the researcher, and the effect on the dependent variable is observed across the resulting time periods graphically through `dynardl.simulation.plot()`.

This functionality simplifies the production and interpretation of dynamic models to a broader class of users. Basic interpretation of dynamic models usually involves calculating the short-run effect of a one-unit increase in an independent variable, as well as the long-run effect. The former is easy to estimate and observe (represented by a $\hat\beta$ coefficient in an ARDL model), but the latter---and its associated confidence intervals---are harder to obtain, since they involve a non-linear combination of parameter estimates. Inferences can also be obscured due to the complexity of the model (like including multiple lags and/or first-differences of dependent and independent variables). For short-run effects, increased model complexity encourages a focus on "boring" effects (the one-unit effect reported by the $\hat\beta$ in the regression output); for long-run effects, greater model complexity results in increasingly intractable closed-form estimates of the effects. Stochastic simulations can solve both of these problems by allowing for short-run dynamics beyond the traditional one-unit effect, and also by removing the need for analytic calculations when discussing long-run effects. 

# Plotting functionality
dynamac implements six new visualizations. These can be presented one at a time using the command `dynardl.simulation.plot()`, or all six using `dynardl.all.plots()`. The resulting plots are:

* The response path of the level of the dependent variable, given a change (a "shock", such as a one standard deviation increase or decrease) in a single independent variable at one point in time.
* Deviations between the predicted and average values of $y_t$: $\hat{y}_t - \bar{y}$. While this plot follows a response path identical to the plot above, if a shock dissipates and the series reverts to its mean, this is easier to observe if we do not have to subtract the stable starting value from the plot. If the series reverts to a value other than its mean, we would be especially interested in this new long-run mean in response to the shock. This is analogous to an impulse response function.
* The period-to-period response in $y_t$ by differencing the response path. This allows us to visualize successive movements over time.
* The size of the shock itself over time. This allows us to visualize dynamic persistence, so we can make statements like "how many time points until the shock is at 50 percent of its original magnitude?" or "how many time points until the shock is effectively zero?" This information is poorly represented by either the short-run or long-run effects as traditionally calculated.
* The cumulative response in $y_t$ to a shock in $x_t$. This helps visualize what the "final" effect of a shock is on the dependent variable. A variable might first cause a positive response, but ultimately may be overwhelmed by a negative movement. Unlike other plots, this considers the cumulative histories in each individual simulation when plotting the response in $y$ over time, making the confidence intervals more conservative. This is analogous to the traditional long-run effect of an ARDL model.
* Total movement that occurs in $y_t$ given a change in $x_t$, which is analogous to the "absolute" effect of a shock to $x_t$ on $y_t$. In this setup, each of the period-over-period changes sums towards the "total" movement in $y_t$, so that all of the movement in $y$ is attributed as an effect to the independent variable.

![Six quantities of interest from the ARDL equation $\Delta y_t = -0.8 y_{t-1} -2 \Delta x_t  + x_{t-1} + u_t$.](allplots-revised.pdf)

Figure 1 plots these visualizations from a hypothetical ARDL model. Starting with the top-left plot and moving clockwise, users can easily understand how a dependent variable responds over time to a shock in an independent variable, how it changes away from its average pre-shock value, how the change itself changes over time, the cumulative absolute changes in the dependent variable, the cumulative nature of the changes in the dependent variable, and how the shock decays over time. All six graphical quantities of interest provide a richer set of information about our inferences. These require no additional mathematical burden on the user, regardless of the underlying model complexity. 

# Conclusion

Dynamic models are complicated since they have effects that do not lend themselves to straightforward interpretation. Stochastic simulations can move our interpretations forward, even in the face of increased model complexity. As we have shown, there are six quantities of interest that are of interest to users, since they provide us with a richer understanding of the dynamic process taking place. dynamac now includes the ability to obtain all of these visualizations described above in a single function: ``dynardl.all.plots()``. This will help users expand the different types of quantities of interest they can use to assess statistical and substantive significance.

# References

