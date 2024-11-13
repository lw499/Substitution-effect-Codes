This repository contains code from the paper "Identification and estimation of the average causal effects under dietary substitution strategies", co-authored by Chiu and Wen. 
The files contain codes to generate data used in simulation studies and methods to estimate causal estimands described in the main paper. 
The files to reproduce the results found in the main paper (and in the supplementary materials) include

1. datagen.R: contains the function to generate data sets
2. deterministic_datagen_wide_true.R: code to produce the true parameter estimates
3. dr.R: codes to produce the parameter estimates and standard errors from the doubly robust method described in the paper
4. ipw.R/gcomp.R: codes to produce the parameter estimates and standard errors from the singly robust estimators described in the paper
