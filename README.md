## Description

This repository accompanies the manuscript _Bayesian Inference of Mixing and Transmission Heterogeneity in Stratified Disease Surveillance Models_. Section 3 of the manuscript details a simulation study whose code is included in the `simulation-study` folder, and Section 4 details a model application whose code is included in the `noro-inference` folder. Code provided for the COVID-19 application described in Section 5 is intentionally left incomplete: in-accordance with our data use agreement, reproduction of these results is not available at this time.

## Publication

ArXiV: 

## Tips & Tricks

1. If running this on your computing cluster, consider making a `model.stan` and a separate `model-nogq.stan`, where the latter has no generated quantities block. The GQs for these models are very expensive to store. With the two files separate, you can run `model-gq.stan` first and run the generated quantities block of `model.stan` separately-- this is why [script.R] is written to delete GQ draws immediately after summarizing them. 
