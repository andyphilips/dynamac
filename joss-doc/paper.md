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
date: 9 April 2020 
tags: 
  - time series
  - methodology
  - long-run effects
  - simulations
  - quantities of interest
bibliography: jordan-philips.bib 
---

# Summary 

dynamac [-@dynamacCRAN] uses stochastic simulation techniques [@jordan2018cointegration;@jordan2018dynamic] to easily recover traditional quantities of interest, such as short- and long-run effects, even from complex dynamic specifications (such as the addition of multiple lags of dependent variables, or independent variables appearing in multiple lags or multiple forms like contemporaneous and lagged differences, or contemporaneous and lagged levels), including models with cointegration [@philips2016have]. dynamac implements these tools in a cohesive ecosystem of estimating and interpreting autoregressive distributed lag (ARDL) models. The procedure is as follows. First, the parameters are estimated using an ARDL model. Next, the parameters used to draw $s$ simulations from a multivariate normal distribution, with mean $\hat{\beta}$ and variance obtained from the estimated variance-covariance matrix. Stable predicted values for $y$ using the $s$ simulations of $\hat{\beta}$ are then created, using starting values of the independent variables (usually means for continuous variables or modes for categorical ones). Users can choose from either expected values or---by incorporating fundamental model uncertainty---predicted values. Next, at a defined time period $t$, one of the independent variables is "shocked" by an amount determined by the researcher, and the effect on the dependent variable is observed across the resulting time periods. To summarize, models are estimated and counterfactual shocks are simulated through `dynardl()`. Then post-estimation dynamic inferences can be uncovered through the new command `dynardl.simulation.plot()`, which has six plotting responses mapped to the six visualizations that we outlined below. 

This functionality vastly simplifies the production and interpretation of dynamic models, increasing their applicability to a broader class of users. Dynamic models, like the ARDL model, can provide researchers with a variety of inferences about both short- and long-run effects of variables. Basic interpretation of dynamic models usually defaults to presenting both the short-run effect of a one-unit increase in an independent variable on the dependent variable, as well as the long-run, or cumulative total effect. The short-run effect is relatively easy to estimate and observe (represented by a $\hat\beta$ coefficient in an ARDL model), but long-run effects---and their associated confidence intervals---are harder to obtain, since they involve a non-linear combination of parameter estimates.  As a result, they are often ignored, even though we are meant to be taking them "seriously" [@de2008taking].

In addition to the difficulties in calculating short- and long-run effects, inferences can also be obscured due to the complexity of the model. This is especially true as the models become more sophisticated (like including multiple lags and first-differences of dependent and independent variables). For short-run effects, increased model complexity encourages a focus on "boring" effects (the one-unit effect reported by the $\hat\beta$ in the regression output); for long-run effects, greater model complexity results in increasingly intractable closed-form estimates of the effects. Stochastic simulations can solve both of these problems by allowing for short-run dynamics beyond the traditional one-unit effect, and also by removing the need for analytic calculations when discussing long-run effects. 

# Plotting functionality
Specifically, dynamac implements six new visualizations. These can be presented individually using `dynardl.simulation.plot()`, or all together using the command `dynardl.all.plots()`:

* The response path of the level of the dependent variable, given a change (a "shock") in a single independent variable at a point in time [@williams2012but]. Note the value of this "shock" is user-defined, allowing for meaningful short-run inferences beyond the typical one-unit effect represented by a $\hat\beta$; for instance, a standard deviation increase or decrease.
* Deviations between the predicted and average values of $y_t$ (i.e., $\hat{y}_t - \bar{y}$). This is analogous to an impulse response function. If a shock dissipates and the series reverts to its mean, this is easier to observe if we do not have to subtract the stable starting value from the plot. If the series reverts to a value other than its mean, we would be especially interested in this new long-run mean in response to the shock.
* The period-to-period response in $y_t$ by differencing the response path. Such a figure allows us to visualize successive movements over time; if the change is equal to zero, it means there is no movement between the past and current period.
* The cumulative response in $y_t$ to a shock in $x_t$. This helps us to understand what the "final" effect of a shock to an independent variable does to a dependent variable. An independent variable might first cause a positive response, but ultimately may be overwhelmed by a negative movement as time goes on. This is analogous to the traditional long-run effect of an ARDL model; however, it requires no analytical calculation.
* Total movement that occurs in $y_t$ given a change in $x_t$. We might think of the "absolute" effect of a shock to $x_t$ on $y_t$ as being the sum over these positive and negative effects, but without allowing them to offset. In this setup, each of the period-over-period changes counts towards the "total" movement in $y_t$, so that all of the movement in $y$ is attributed as an effect to the independent variable.
* The size of the shock itself over time. This allows us to visualize dynamic persistence, so we can make statements like "how many time points until the shock is at 50 percent of its original magnitude?" or "how many time points until the shock is effectively zero?" This information is poorly represented by either the short-run or long-run effects as traditionally calculated, but is of interest to analysts.

![Six quantities of interest from the ARDL equation $\Delta y_t = -0.8 y_{t-1} -2 \Delta x_t  + x_{t-1} + u_t$.](allplots.pdf)

Figure 1 plots these visualizations from a hypothetical ARDL model. Users can easily understand how a dependent variable responds over time to a shock in an independent variable, how it changes away from its average pre-shock value, how that change itself changes over time, how the shock decays over time, the cumulative nature of the changes in the dependent variable, and the absolute changes in the dependent variable. All six graphical quantities of interest provide a richer set of information about our inferences.

The additional nicety is that this requires no additional mathematical burden on the user, regardless of the underlying complexity of the model. This is critical, as we're often encouraged to make modeling decisions that increase complexity, not realizing the resulting strain on analytic inferences. For instance, Wilkins [-@wilkins2018lag] advocates adding an additional lag of $y_t$ at $t-2$ in order to remove residual autocorrelation. Following this advice, the long-run effect is more difficult to calculate, as it is now a non-linear function of both coefficients on $x_t$ at $t$ and $t-1$ as well as the coefficients on the lag of $y_t$ and the additional lag at $t-2$. Such complexity could reasonably encourage a practitioner to ignore long-run inferences; this is obviously not the preferred outcome. Through dynamac, though, stochastic simulations retrieve these quantities with no analytic calculations required. And the resulting graphic inferences are easy to both create and interpret, allowing users to focus on what they know best: their underlying models.

# Conclusion

Dynamic models are complicated since they have many "hidden" effects that do not lend themselves to straightforward interpretation. Scholars that use these models typically interpret only the short- and long-run effects. Stochastic simulations can move our interpretations forward, even in the face of increased model complexity. As we have shown, there are six quantities of interest that are of interest to users, since they provide us with a richer understanding of the dynamic process taking place. dynamac now includes the ability to obtain all six of the visualizations above in a single function: ``dynardl.all.plots()``. In doing so, we help users expand the different types of quantities of interest they can use to assess statistical and substantive significance of their independent variables. Such visualizations should find broad usefulness, given the usage of dynamac for applied questions in fields from environmental science [@sotte2019;@sotte2020;@jcp,@ijerph] to economics [@espr;@te] to political science [@trumpy] and beyond.

# References
