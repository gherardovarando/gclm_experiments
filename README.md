# Simulations and Experiments for Graphical Continuous Lyapunov Models (2020)  

This repository contains the files for replicating the simulations
and the examples in Graphical Continuous Lyapunov Models, G. Varando and N. R. Hansen 
2020.   

## Required packages

* `glmnet`
* `glasso`
* `igraph`
* `MASS`
* `ggplot2`

## Simulations

### Run the simulations

To run the simulations just move in the directory and run the R script 
`simulate.R`. 
The simulations can take up to a couple of days, depending on the
machine.
The exact simulations run in the paper are the ones 
that can be executed from the script `run_simulate.sh`. 

### Plot results 

Before plotting the results you need to compute the evaluations and 
aggregate the results with the R script `aggregate_results.R` or directly 
using the script `run_aggreagate.sh`.  

Then you can lunch the plotting scripts `plot_simulated.R`,
`plot_marginalized.R` and `plot_times.R` 
 to generate the same plots as in the paper.  

## Protein-signaling network

First run the script `proteins.R` to generate the results. 

The estimated graphs can be obtained running `plot_proteins.R`. 

