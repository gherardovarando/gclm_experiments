# Simulations and Experiments for TITLE HERE 

This repository contains the files for replicating the simulations
and the example in TITLE 

## Required packages

* `glmnet`
* `igraph`
* `MASS`
* `ggplot2`

## Simulations

### Run the simulations

To run the simulations just move in the directory and run the R script 
`simulate_ou.R`. 
The simulations can take up to a couple of days, depending on the
machine.
If you have access to a computing machine you can lunch parallel jobs 
using the script `run_simulate.sh`.  

### Plot results 

Before plotting the results you need to compute the evaluations and 
aggregate the results with the R script `aggregate_results.R`. 

Then you can lunch the plotting scripts `plot_simulated_ou.R`
and `plot_avgprc.R` to generate the same plots in the paper.  

## Protein-signaling network

First run the script `proteins.R` (or lunch parallele computations with 
`run_proteins.sh`) to generate the results. 

The estimated graph can be optained running `plot_proteins.R`. 

