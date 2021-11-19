logisticLAGO Package
================

## Description

The “Learn-As-You-Go” (LAGO) design is motivated by large scale public
health intervention trials, where an intervention package is
prospectively modified over each stage such that an “*optimal*”
intervention package is rolled out to the participants in the next
stage. This R package implements the optimization for the LAGO study
with a single constraint, namely the outcome goal constraint and binary
response.

## Installation

You can install logisticLAGO R package from GitHub with:

``` r
 devtools::install_github("Arhit-Chakrabarti/logisticLAGO")
```

Note that the library *devtools* needs to be installed before installing
the R package from GitHub.

## Usage

Once the **logisticLAGO** library is installed load the library in the R
workspace.

``` r
library("logisticLAGO")
```

The basic setup of the LAGO design requires a multi-component,
continuous intervention package and a binary response variable which may
be for example, how well the participants perform on a test following
administration of the intervention package. The intervention package is
prospectively changed at every stage based on the cumulative data
collected over the stages, such that an *optimal* package is rolled out
to the participants in the next stage such that the probability of
success for the binary response is above a desired threshold while
minimizing the implementation costs.

### In a ongoing trial

In the case when a trial has already been designed and data from the
trial has been collected, the next step of the study is to estimate the
optimal intervention package to be rolled out to the participants in the
next stage such that the otcome goal of the study is met and costs
minimized. To estimate the optimal intervention package the vector of
per unit linear costs for the intervention package (*cost\_lin*), the
vector of minimum (*x.l*) and maximum (*x.u*) values of the components
of the intervention package, desired outcome goal (*p\_bar*), the
estimated *β̂* from fitting a logistic regression model to the observed
response, which gives the estimated effect of the corresponding
intervention package and the intervention package rolled out at the
current stage (*x.init*). The estimated optimal interventio is given by

``` r
opt_lago = opt_int(cost = cost_lin, beta = beta, lower = x.l, upper = x.u, pstar = p_bar, starting.value = x.init)
```

### Before starting a trial

#### Single center LAGO design

This package may also be used before starting a trial to get an idea
about the optimal intervention package based on an initial intervention
package and best guesses about the effects of the components on the
response. Initial package and idea about the effects of the components
of the intervention package may be obtained from an investigator or from
knowledge of prior or concurrent intervention trials. The simplest case
is when the trial is designed to be conducted in a single center or
location. The number of stages (*K*) in the LAGO design, sample size per
stage (*n*), the unit costs for the intervention package components
(*cost\_lin*), the vector of minimum (*x.l*) and maximum (*x.u*) values
of the components of the intervention package, desired outcome goal
(*p\_bar*), the best guess intervention effect (*beta*) and the initial
value of intervention package (*x.init*) along with the expected
variation in rolling out the intervention package (*icc*) are used to
simulate the estimated intervention package over the different stages,
the estimated outcome goal and the estimated power of the test of
“*no-intervention*” effect at the end of the study. This is done using
the following function:

``` r
sim_sc <- sc_lago(x0 = x.init, lower = x.l, upper = x.u, nstages = K, beta.true = beta, sample.size = n, icc = icc, cost.vec = cost_lin, prob = p_bar, B = 100, intercept = TRUE)
```

#### Multi-center LAGO design

Another common study design is when the trial is planned to be conducted
in a multiple centers or location. Apart from the information as
required by the single center design, the number of centers per stage
(*J*) and sample size per center per stage (*njk*) is used to simulate
the optimal intervention package to be rolled out in each of the centers
in the next stage such that the goals of the study are met. As before,
estimated power for the test of “*no-intervention*” effect at the end of
the study is also provided through simulations. The corresponding
function is:

``` r
sim_mc <- mc_lago(x0 = x.init, lower = x.l, upper = x.u, beta.true = beta, nstages = K, centers = J, sample.size = njk, icc = 0.1, prob = p_bar, cost.vec = cost_lin)
```

## Details

For more information on logisticLAGO Package, please access the package
documentations or vignettes. Please feel free to contact the author.
