# GMMuni
This repository contains the R functions and the replicability file used for the paper "Unimodal density-based clustering and merging algorithm using Gaussian mixtures" by

D. Tancini (University of Perugia, IT);
L. Scrucca (University of Bologna, IT);
F. Bartolucci (University of Perugia, IT).

# ----------------------------
Replicability file: Replicability.R
# ----------------------------

This file contains:

- Motivating example 1;
- Motivating example 2;
- Motivating example 3;
- A subset of the simulation design for the merging problem;
- A subset of the simulation design for the unimodal problem;
- Simulation design when non informative features are considered;
- Application Bankruptcy;
- Application Thyroid.

# ----------------------------
Function Name: GMMuni
# ----------------------------

Description: The GMMuni function is designed to estimate unimodal GMMs.

Input: 

- x, a dataset of dimension n x d, where n stands for the number of observations and d as the dimension;

- R, the max number of iterations of the algorithm;

- G, the number of components;

- gamma, a tuning parameter for the penalization function (see 4.2 in GMM_Unimod.pdf Dropbox);

- set, if 1 => m is defined using the mode of the initial GMM estimated, else => m is defined using the weights of the initial GMM estimated 

- model, if "yes" a Mclust model is used for evaluate starting points and mode estimates

- model_est, Mclust model

- force_est, if 1 the algorithm is forced to provide a PMLE even if MEM indicates 1 mode

- den, a logical, if TRUE a denoising procedure is used when d>1 to discard all modes whose density is negligible;

Output:

- Mus, the estimated means;

- Sigmas, the estimated variance-covariances;

- pis, the estimated weights;

- classification, the associated component for each observation (MAP);

- m, mode hyperparameter;

- MuMclust, the initial means estimated using Mclust;

- SigmaMclust, the initial variance-covariances estimated using Mclust;

- PiMclust, the initial weights estimated using Mclust;

- p2, the value of the penalization component;

- n_modes, the number of modes estimated using a MEM with Mus, Sigmas, pis;

- x, the initial dataset.

# ----------------------------
Function Name: Hellingher_algorithm
# ----------------------------

Description: The Hellingher_algorithm function is designed to merge Gaussian components.

Input: 

- mod, a Mclust model estimated 

- ep, \epsilon value defined in Hellinger_mixture_algorithm

Output:

- merg, the components to be merged

- H, the estimated Hellinger matrix 
