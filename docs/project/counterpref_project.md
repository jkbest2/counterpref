---
title: "Designing a survey around fishery-dependent data"
subtitle: "FISH 507B final project"
author: John Best
bibliography: ../fish507b.bib
header-includes: \renewcommand{\vec}[1]{\boldsymbol{#1}}
---

# Introduction

Model-based indices of abundance are increasing in popularity, particularly
those that use spatiotemporal statistical models for standardization
[@ThorsonEtAl2020]. When available, these spatiotemporal model-based indices are
usually based on fishery-independent observations, even when there are a
substantial number of fishery-dependent observations available. This fails to
utilize one of the advantages of model-based indices: the ability to incorporate
multiple sources of data. Many fisheries have extensive observer programs that
record detailed fishery-dependent catch information. Naturally, there are many
more fishery-dependent observations than could feasibly be collected during a
survey. Incorporating these fishery-dependent observations into indices of
abundance has the potential to increase the precision of these indices.

Moving from an index of abundance that uses fishery-independent data exclusively
to one that also includes fishery-dependent data is not without pitfalls.
Fishery-dependent data are not collected using a defined sampling procedure.
This may result in higher variance of catches (e.g. among vessels). Fishing
locations are also chosen for their potential profit rather according to a
specified design. This targeting behavior may result in a preferential sample,
where only areas with higher biomass of the targeted species are sampled by the
fishery. This may bias resulting estimates of abundance [@GelfandEtAl2012;
@Conn2017]. It is possible to account for the effects of preferential sampling
post hoc, but this comes at the cost of considerable model complexity
[@Diggle2010; @DinsdaleSalibian-Barrera2019a]. A simpler option (from a modeling
perspective though admittedly not a logistical perspective) would be to design a
survey that accounts for the weaknesses of the fishery-dependent observations.

At one extreme, survey locations could be chosen so that overall sampling effort
is roughly uniform throughout the domain, eliminating the observation location
preference. This would generally require an excessive amount of survey effort. A
second option would be to overlay a standard survey grid over the domain of
interest, as is currently done. This would ensure samples outside the area where
fishery-dependent data is abundant, but apportions survey effort in a way that
does not account for the fishery-dependent samples. In between these options, a
survey design with greater sampling intensity in under-fished areas can fill in
information where it is most needed. Any successful survey that intends to
integrate fishery-dependent and -independent observations will need to allow
enough spatial overlap between the two so that differences in catchability may
be estimated. Sampling away from a species' area of highest abundance can also
aid in establishing the limits of a species' range, and provide additional
information about species that are not the target of commercial fisheries.
Survey design under preferential sampling has been studied in the case where the
preference is meant to be preserved [@daSilvaFerreiraGamerman2015], designing a
survey specifically to *counteract* the effects of preferential sampling has not
been addressed.

There are two key questions. The first is how the number of
counter-preferentially chosen survey stations affects the bias of an index of
abundance estimated with varying numbers of fishery-dependent observations. The
second question is how well differences in catchability between the
fishery-dependent and -independent vessels can be estimated when the survey
stations are chosen counter-preferentially and thus will have less overlap with
the fishery-dependent data.

# Methods

## Operating model 

A single simulation study will be used to evaluate the two questions presented
above. Catches were be simulated using the
[`FisherySim.jl`](https://github.com/jkbest2/FisherySim.jl) software package.
The simulated fishery occurs on a $100\times100$ grid, where each cell is
assigned a biomass. Total initial biomass is fixed; the spatial distribution of
this biomass is determined by a spatially correlated habitat covariate and a
habitat preference function. This habitat preference, along with a distance
penalty, determine a movement operator. The initial spatial distribution of the
population is the equilibrium distribution under this movement operator (i.e.
the eigenvector associated with its largest eigenvalue). Figure
{@fig:init_abund} is an example of an initial abundance distribution.

![One replicate of an initial abundance distribution.](../../figs/initial_abundance.svg){#fig:init_abund}

Given the initial abundance distribution, the first step in the simulation is
biomass removal by fishing. The fishery-dependent fleet targets 2000 cells (with
replacement), where the probability of choosing cell $s$, $p_{comm}(s)$ is the
ratio of the biomass in that cell, $n(s)$, to the total biomass

$$p_{comm}(s) = \frac{n(s)}{\sum_s n(s)}.$$

Fishery-dependent fishing locations are chosen at random at the start of each
simulation year. Figure {@fig:comm_effort} is an example of the spatial
allocation of commercial fishing effort given the abundance distribution
illustrated in figure {@fig:init_abund}.

![Example of commercial fishing effort under a linear preference.](../../figs/commercial_effort.svg){#fig:comm_effort}

Fishery-independent survey stations were also chosen at random (without
replacement) from a 400-station regular grid over the spatial domain. Depending
on the scenario, 50, 100, 200, 300, or all 400 stations were chosen. The
probaiblity $p_{surv}(s)$ of a station being chosen was based on the sum of the
abundance of the 25 nearest cells, $n^*(s)$. Weights were inverted so that
stations with lower abundance had a higher probability of being chosen, based on
the overall mean $\bar{n}^*$, plus a constant $k$ so that all of the weights are
positive[^1].

$$p_{surv}(s) = 2 \bar{n}^* - n^*(s) + k$$

[^1]: Future iterations of this study will use the [softmax
    function](https://en.wikipedia.org/wiki/Softmax_function), so that the
    counter-preferential weights can just be the negated abundance values, so
    that $$p_{comm}(s) = \frac{\exp [n(s)]}{\sum_s \exp [n(s)]}$$ and
    $$p_{surv}(s) = \frac{\exp [-n^*(s)]}{\sum_s \exp [-n^*(s)]}.$$
    
Figure {@fig:surv_stations} shows the probability that each survey station will
be chosen given the initial population distribution in figure {@fig:init_abund}.
These stations are chosen at the beginning of the simulation and used for each
simulated year.

![Probability that a survey station will be used given the abundance
distribution in figure {@fig:init_abund}](../../figs/station_prob.svg){#fig:surv_stations}

Fishery-dependent fishing events occur in random order, while the
fishery-independent events occur in a fixed order, from the lower left to upper
right corner of the spatial domain. Fishery-dependent and -independent fishing
events occur in random order. Catchability for both vessels is fixed at 0.2, and
cells may become locally depleted within each season. Catch observations are
simulated using a compound Poisson gamma distribution, a special case of the
Tweedie distribution [@Tweedie1984; @Shono2008]. This allows for catches that
are exactly zero or positive. In this case the Tweedi dispersion parameter
$\phi$ is fixed to $1.2$ and the Tweedie shape parameter $p$ is fixed to $1.84$.

After harvest, the population produces additional biomass according to a
Beverton-Holt surplus production function. Carrying capacity is region-wide, and
fixed at 100 for all simulations. Produced biomass is distributed
proportional to the current abundance. For total abundance $N_t = \sum_s n_t(s)$,

$$n_{t+1}(s) = \frac{n_t(s)}{N_t} r N_t \left(1 - \frac{N_t}{K}\right).$$

Finally, biomass is redistributed via the movement operator described above,
based on habitat preference and a distance penalty. Once this is complete, a new
year begins and harvest can occur. True abundance indices are calculated based
on the abundance at the start of each year. 25 habitats were generated, and then
fisheries were simulated for 15 years with 50, 100, 200, 300, or 400 survey stations.
    
## Estimation model

Indices of abundance were estimated using the
[`spatq`](https://github.com/jkbest2/spatq) R package, which includes a
spatiotemporal index standardization model written in Template Model Builder
[@KristensenEtAl2016]. A matching compound Poisson gamma observation likelihood
is used, where for catch $y_i$, expected catch $c_t(s)$ at time $t$ and location $s$,

$$y_i \sim \operatorname{Tweedie}(c_t(s), \phi, p).$$

Expected catch is modeled using the linear predictor

$$\log c_t(s) = \vec{X}\vec{\beta} + \vec{A}_\omega \vec{\omega} + \vec{\Gamma}\vec{\lambda},$$

where $\vec{X}$ is the design matrix for the year effect, $\vec{\beta}$ is the
vector of year effects, $\vec{A}$ is the spatial projection matrix,
$\vec{\omega}$ is the vector of spatial random effects, and $\vec{\Gamma}$ is
the design matrix of the vessel effect and $\vec{\lambda}$ is the vessel
effect. The spatial random effect takes the distribution

$$\vec{\omega} \sim \operatorname{MVN}(\vec{0}, \vec{Q}^{-1}).$$

The precision matrix $\vec{Q}$ is constructed to approximate a Mat\'{e}rn
Gaussian random field with smoothness 1, using the SPDE approximation
[@Lindgren2011].

The index of abundance for each year, $\hat{I}_y$ is estimated by using the
design matrix for year $y$ $\vec{X}^*_y$, the estimated year effects
$\hat{\vec{\beta}}$, and projecting the spatial random effect to a regular grid
using the matrix $\vec{A}^*_\omega$ so that the abundance can be summed across
the spatial domain. In this case the integration weights can be set to one:

$$\hat{I}_y = \vec{1}^{T} \exp(\vec{X}^*_y \hat{\vec{\beta}} + \vec{A}^*_\omega \hat{\vec{\omega}}).$$

Standard errors for the index estimates are estimated using the delta method
with bias correction [@ThorsonKristensen2016]. A total of 15 estimation models
were fit for each of the 25 replicates. These used each combination of 50, 100,
200, 300, or 400 survey stations and 500, 1000, or 2000 fishery-dependent
observations (selected at random from the 2000 total observations).

## Evaluation

Bias is an important component of the error of an index of abundance. a log
regression of the estimated index on the true total abundance $b_{roey}$ of
replicate $r$ with operating model $o$, estimation model $e$, andy year $y$,
with appropriate offsets provides a bias estimate $\delta_{oe}$ [@Thorson2015]:

$$\begin{aligned}
\log \hat{I}_{roey} &= \alpha_{roe} + \delta_{oe} \log b_{roey} +
\varepsilon_{roev}\\
\varepsilon_{roev} &\sim \operatorname{Normal}(0, \sigma^2_oe).
\end{aligned}$$

The bias metric $\delta$ will be one if the estimates are unbiased.

Total root mean square error of the index estimate can be evaluated after
rescaling the indices and true abundances to a common scale, so that

$$\begin{aligned}
\hat{I}'_t &= \frac{\hat{I}_t}{\exp \sum_t \log \hat{I}_t},\ \text{and}\\
{b}'_t &= \frac{{b}_t}{\exp \sum_t \log {b}_t}.
\end{aligned}$$

This also allows for evaluation of the coverage of the confindence intervals
based on the estimated standard errors of the indices.

Bias and error of the vessel effect are more straightforward to evaluate because
they are naturally on the same scale.

# Results

Note that the following results are based on 25 replicates, and additional
replicates would smooth out some of the noise that may be present.

Bias due to preferential sampling is apparent, particularly when there are many
fishery-dependent observations and few fishery-independent observations. This
bias appears to decrease as additional survey stations are added, but does not
entirely disappear when there are substantially more fishery-dependent
observations (figure @fig:bias_plot2).

![Bias metric for each combination of number of survey stations and commercial
observations, using 15 years from 25 replicates.](../../figs/bias_plot2.svg){#fig:bias_plot2}

Unsurprisingly, root mean square error decreases with additional additional
observations (figure @fig:rmse_plot).

![Root mean square error for each combination of number of survey stations and
commercial observations, using 15 years from 25 replicates.](../../figs/rmse_plot.svg){#fig:rmse_plot}

When fishery-dependent effort is distributed preferentially and survey stations
are allocated counter-preferentially, there will naturally be limited overlap of
observations. This makes estimating catchability differences between the vessel
types difficult. In these simulations, there was no difference in catchability
between the fishery-independent and -dependent vessels. The catchability offset for the fishery-dependent vessels
($\lambda$) should then be zero. There is a clear bias in the estimate of this
parameter when few survey stations are used, but it disappears when 300 and 400
survey stations are used (figure @fig:lambda_val).

![Density plot of catchability offset for the fishery-dependent vessel over 25 replicates. The
generative value, $0$, is marked with the dashed line.](../../figs/lambda_val.svg){#fig:lambda_val}

The standard error of the catchability offset parameter $\lambda$ naturally
decreased with additional observations. Additional survey stations reduced the
standard error more than additional fishery-dependent observations (figure @fig:lambda_sd)

![Density plot of catchability offset standard error estimates over 25 replicates.](../../figs/lambda_sd.svg){#fig:lambda_sd}

# Discussion



## Extensions

To ensure that this project can feasibly be completed in the required timeframe,
a number of interesting extensions will not be addressed. The first issue is
that the focus here is on estimating an index of abundance; the sampling
strategies and evaluation criteria do not account for other data types that may
not be available from fishery-dependent data. Second, practical implementation
of these survey designs would require strategies for specifying survey locations
*before* fishery-dependent effort distribution is known. This may be trivial in
fisheries that don't change much year-to-year, but it may be necessary to
include uncertainty in fishery-dependent location preference in order to make
the strategy robust to spatial shifts in the fishery. Third, this simulation
study focuses on a single species, which exclusively drives the
fishery-dependent fishing location preference. Survey designs in this case may
or may not be optimal when the location preference is driven by other factors
such as a different species. It would also be interesting to see how these
strategies fare when the goal is to produce indices of abundance for multiple
species.

- marginal information for additional preferential commercial observation vs.
  additional survey station
# References
