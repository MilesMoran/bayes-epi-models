## Overview

This repository accompanies the manuscript _Bayesian Inference of Mixing and Transmission Heterogeneity in Stratified Disease Surveillance Models_. Section 3 of the manuscript details a simulation study whose code is included in the `simstudy` folder, and Section 4 details a model application whose code is included in the `inf-noro` folder. Code provided for the COVID-19 application described in Section 5 is intentionally left incomplete under `inf-covid`: in-accordance with our data use agreement, reproduction of these results is not available at this time.

## Detailed Description

Section 3 of the manuscript describes a simulation study comparing the negative-binomial version of our proposed model (see Section 4) to an existing competitor (see Meyer and Held, 2017). We refer to our proposed model as the "full" model and the Meyer-Held comparator as the "reduced" model, in the sense that the Meyer-Held model is retrievable as a special/limiting case of our proposal. 
???

Section 4 of the manuscript describes a negative-binomial version of our proposed model class, in which
(1) the negative-binomial conditional likelihood for incident counts is _not_ adjusted for small population counts (susceptible pool size is always taken to be the population count); 
(2) the mean of the latent-infectiousness variable `r` is scaled ONLY by lag-1 incidence under the rare-disease assumption that prevalence equals incidence; and,
(3) the age-group mixing weights are calculated as the column-normalization of an estimated contact-rate matrix $C$, with prior described in Section 2.
The code used to fit this model has been provided in the `inf-noro` subdirectory. All data are accessed through the `hhh4counts` R package,
????

Section 5 of the manuscript describes a beta-binomial version of our proposed model class, in which 
(1) the beta-binomial conditional likelihood for incident counts is bounded-above by an time-varying estimate of the susceptible pool size (`x_hat` in Stan);
(2) the mean of the latent-infectiousness variable `r` is scaled by a disease prevalence estimate $\hat{Y}_{tgi}$ more complex than the rare-disease assumption that prevalence equals incidence; and,
(3) the age-group mixing weights are calculated as the column-normalization of an estimated contact-rate matrix $C$, with prior described in Section 2.
The code used to fit this model has been provided in the `inf-covid` subdirectory. Note that the original incidence data is omitted. If you generate data with the same schema, the main script `inf-covid/R/mi-covid-fit.R` should work as-intended. This script references static copies of the population counts matrix $E$ and the distance matrix $D$, so the scripts used to compute those are included in `inf-covid/R` and the publicly-available data needed to do so are included in `inf-covid/data`. 

## Publication

ArXiV: 

## Tips & Tricks

1. If running the simulation study code on your computing cluster, consider making a `model.stan` and a separate `model-nogq.stan`, where the latter has no generated quantities block. The GQs for these models are very expensive to store. With the two files separate, you can run `model-gq.stan` first and run the generated quantities block of `model.stan` separately-- this is why [script.R] is written to delete GQ draws immediately after summarizing them. 
