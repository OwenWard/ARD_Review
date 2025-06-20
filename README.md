# Bayesian Modeling for ARD, a Unified Perspective

This repository provides code to simulate 
ARD from two popular ARD models:

1. The Latent Space Model of McCormick and Zheng, 2015

2. The NSUM Model of Maltiel et al, 2015


We provide simulated data from each of these models, along
with Stan code to fit the following popular ARD models.

1. Two null models, considering both fixed and varying degree 
parameters for each node in the network

2. The overdispersed model of Zheng et al, 2006.

3. The Latent Space Model of McCormick and Zheng, 2015.

For each of these models we provide Stan code incorporating 
scaling, similar to that considered by Laga et al (2023).

We also provide R functions to perform posterior predictive checks
on these models.

Finally, we consider an approach to do model selection using K-fold
cross validation. This is not applicable to all potential
ARD models in the literature.

The code to simulate this ARD and fit each of the above models
is included in `Simulation_Scripts`.