---
title: "Designing a survey around fishery-dependent data"
subtitle: "FISH 507B project proposal"
author: John Best
bibliography: fish507b.bib
---

# Introduction

Model-based indices of abundance are increasing in popularity, particularly
those that incorporate spatiotemporal dynamics [@ThorsonEtAl2020]. When
available, model-based indices are usually based on fishery-independent
observations. However, one advantage of model-based indices is the ability to
incorporate multiple sources of data. Many fisheries have extensive observer
programs that record detailed fishery-dependent catch information. Naturally,
there are many more fishery-dependent observations than could feasibly be
collected during a survey. Incorporating these fishery-dependent
observations into indices of abundance has the potential to increase the
precision of these indices.

Moving from a fishery-independent index of abundance to one that also includes
fishery-dependent data is not without pitfalls. Fishery-dependent data are not
collected using a defined sampling procedure. This may result in higher variance
of catches (e.g. among vessels). Fishing locations are also chosen for their
potential profit rather according to a specified design. This targeting behavior
may result in a preferential sample, where only areas with higher biomass of the
targeted species are sampled by the fishery. This may bias resulting estimates
of abundance [@GelfandEtAl2012; @Conn2017]. It is possible to account for the
effects of preferential sampling post hoc, but this comes at the cost of
considerable model complexity [@Diggle2010; @DinsdaleSalibian-Barrera2019a]. A
simpler option (from a modeling perspective though admittedly not a logistical
perspective) would be to design a survey that accounts for the weaknesses of the
fishery-dependent observations.

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

# Methods

A simulation study will be undertaken to compare indices of abundance that
include preferentially sampled fishery dependent data as well as
fishery-independent data sampled according to two different design strategies.
The first, reference strategy, overlays a uniform survey grid over the domain,
similar to current practice. The second chooses survey locations at random,
where locations with higher fishery-dependent fishing effort are less likely to
be chosen for the survey. These will be compared over multiple levels of total
survey effort (e.g. 50, 100, 150, and 200 total trawls).

Catches will be simulated using the
[`FisherySim.jl`](https://github.com/jkbest2/FisherySim.jl) software package.
This package generates fishery observations on a grid. A spatially correlated
habitat covariate determines fish movement. Overall population dynamics follow a
region-wide biomass dynamics model, in this case a Beverton-Holt. Multiple
fleets remove biomass from each cell. Each fleet may have different catchability
and/or targeting behavior. This allows the different survey strategies to be
implemented. To simplify the simulations all of the fishery-dependent fishing
locations will be chosen prior to fishing, rather than applying a dynamic
preference that accounts for within-cell depletion. On a $100\times100$ grid this
should have minimal impact on the resulting catches. Pre-specifying the
fishery-dependent fishing locations allows the survey locations to be
pre-specified as well.

Indices of abundance will be estimated using the
[`spatq`](https://github.com/jkbest2/spatq) R package, which includes a
spatiotemporal index standardization model written in Template Model Builder
[@KristensenEtAl2016]. Indices will be evaluated for bias and total error,
particularly relative to total survey effort.


# Extensions
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

# References
