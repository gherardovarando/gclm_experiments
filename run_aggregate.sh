#!/bin/bash
tmux new-session -d  "Rscript aggregate_results.R path simulations/ p 10 N 1000 P 10 12 15 20 25 30"

tmux new-session -d  "Rscript aggregate_results.R path simulations/ p 20 P 20"
tmux new-session -d  "Rscript aggregate_results.R path simulations/ p 30 P 30"
#tmux new-session -d  "Rscript aggregate_results.R path simulations/ p 40 P 40"
#tmux new-session -d  "Rscript aggregate_results.R path simulations/ p 50 P 50"

#tmux new-session -d  "Rscript aggregate_results.R path simulations/ p 10 N 100 P 10 12 15 20 25 30"

#tmux new-session -d  "Rscript aggregate_results.R path simulations/ p 20 N 100 P 20"
#tmux new-session -d  "Rscript aggregate_results.R path simulations/ p 30 N 100 P 30"
#tmux new-session -d  "Rscript aggregate_results.R path simulations/ p 40 N 100 P 40"
#tmux new-session -d  "Rscript aggregate_results.R path simulations/ p 50 N 100 P 50"
