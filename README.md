logisticLAGO Package
================
Arhit Chakrabarti

-   [Description](#description)
-   [Installation](#installation)
-   [Usage](#usage)
    -   [In an ongoing trial](#in-an-ongoing-trial)
    -   [Before starting a trial](#before-starting-a-trial)
        -   [Single center LAGO design](#single-center-lago-design)
        -   [Multi center LAGO design](#multi-center-lago-design)
            -   [Equal number of centers per stage and equal sample size
                in each
                center](#equal-number-of-centers-per-stage-and-equal-sample-size-in-each-center)
            -   [Unequal number of centers per stage and equal sample
                size in each
                center](#unequal-number-of-centers-per-stage-and-equal-sample-size-in-each-center)
            -   [Unequal number of centers per stage and unequal sample
                size in each
                center](#unequal-number-of-centers-per-stage-and-unequal-sample-size-in-each-center)
-   [Details](#details)

# Description

The “Learn-As-You-Go” (LAGO) design is motivated by large scale public
health intervention trials, where an intervention package is
prospectively modified over each stage such that an “*optimal*”
intervention package is rolled out to the participants in the next
stage. This R package implements the optimization for the LAGO study
with a single constraint, namely the outcome goal constraint and binary
response.

# Installation

You can install logisticLAGO R package from GitHub with:

``` r
# install.packages("devtools")
# install the package without building vignettes
devtools::install_github("Arhit-Chakrabarti/logisticLAGO")
```

**Note**: The library *devtools* needs to be installed before installing
the R package from GitHub. Please note that this package includes
vignettes which provides further guidelines on how to use the functions
in this R package. To install the logisticLAGO package with the
vignettes, use the code

``` r
# install.packages("devtools")
# install the package building vignettes
devtools::install_github("Arhit-Chakrabarti/logisticLAGO", build_vignettes = TRUE, force = TRUE)
```

# Usage

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
collected over the stages and an *optimal* package is rolled out to the
participants in the next stage such that the probability of success for
the binary response is above a desired threshold, while simultaneously
minimizing the implementation costs.

## In an ongoing trial

In the case when a trial has already been designed and data from the
trial has been collected, the next step of the study is to estimate the
optimum intervention package to be rolled out to the participants in the
next stage such that the outcome goal of the study is met and costs
minimized. To estimate the optimum intervention package, the vector of
per unit costs for the components of the intervention package (*cost*),
the vector of minimum (*lower*) and maximum (*upper*) values of the
components of the intervention package, desired outcome goal (*pstar*),
the estimated (*beta*) from fitting a logistic regression model to the
observed response, which gives the estimated effect of the corresponding
intervention package and the intervention package rolled out at the
current stage (*starting.value*) are required. The estimated optimal
intervention is given by the function *opt\_int*

``` r
opt_int(cost = cost, beta = beta, lower = lower, upper = upper, pstar = pstar, starting.value = starting.value)
```

## Before starting a trial

### Single center LAGO design

This package may also be used before starting a trial to get an idea
about the optimal intervention package based on an initial intervention
package and best guesses about the effects of the components on the
response. Initial package and idea about the effects of the components
of the intervention package may be obtained from an investigator or from
knowledge of prior or concurrent intervention trials. The simplest case
is when the trial is designed to be conducted in a single center or
location. The number of stages (*nstages*) in the LAGO design, sample
size per stage (*sample.size*), the unit costs for the intervention
package components (*cost.vec*), the vector of minimum (*lower*) and
maximum (*upper*) values of the components of the intervention package,
desired outcome goal (*prob*), the best guess intervention effect
(*beta.true*) and the initial value of intervention package (*x0*) along
with the expected variation in rolling out the intervention package
(*icc*) are used to simulate the estimated intervention package over the
different stages, the estimated outcome goal and the estimated power of
the test of “*no-intervention*” effect at the end of the study. This is
done using the following function:

``` r
sc_lago(x0 = x0, lower = lower, upper = upper, nstages = nstages, beta.true = beta.true, sample.size = sample.size, icc = icc, cost.vec = cost.vec, prob = prob)
```

### Multi center LAGO design

#### Equal number of centers per stage and equal sample size in each center

Another common study design is when the trial is planned to be conducted
in a multiple centers or location. Apart from the information as
required by the single center design, the number of centers per stage
(*centers*) and sample size per center per stage (*sample.size*) is used
to simulate the optimal intervention package to be rolled out in each of
the centers in the next stage such that the goals of the study are met.
As before, estimated power for the test of “*no-intervention*” effect at
the end of the study is also provided through simulations. The simplest
design in this case, assumes an equal number of centers per stage and
equal number of samples per center in each stage. The corresponding
function is:

``` r
mc_lago(x0 = x0, lower = lower, upper = upper, beta.true = beta.true, nstages = nstages, centers = centers, sample.size = sample.size, cost.vec = cost.vec, prob = prob, icc = icc)
```

#### Unequal number of centers per stage and equal sample size in each center

Another common design under multi-center LAGO study is considering
different number of centers in each stage. However every center at any
stage have equal samples. The corresponding function is:

``` r
mc_lago_uc(x0 = x0, lower = lower, upper = upper, beta.true = beta.true, nstages = nstages, centers = centers, sample.size = sample.size, cost.vec = cost.vec, prob = prob, icc = icc, bcc = bcc)
```

Here the function argument *icc* denotes the expected variation among
subjects within each center at any stage and the argument *bcc* denotes
the expected variation between the centers at any stage while
implementing the intervention package. Thus, this function separates the
between center and within center variability, while generating the data
for simulations.

#### Unequal number of centers per stage and unequal sample size in each center

The traditional design under multi-center LAGO study is considering
different number of centers in each stage as well as different sample
sizes in each center of any given stage. This design arises quite
naturally, as it is intuitive to have both small number of centers and
small sample size at each center at the beginning of the study and
progressively increase the number of centers and sample size as the
study progresses. The corresponding function is:

``` r
mc_lago_uc.us(x0 = x0, lower = lower, upper = upper, beta.true = beta.true, nstages = nstages, centers = centers, sample.size = sample.size, cost.vec = cost.vec, prob = prob, icc = icc, bcc = bcc)
```

Here, the argument *sample.size* is a list denoting the sample size in
each center of different stages.

# Details

For more information on logisticLAGO Package, please access the package
documentations or vignettes. Please feel free to contact the author.
