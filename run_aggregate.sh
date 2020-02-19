#!/bin/bash
tmux new-session -d  "Rscript aggregate_results.R path simulations0/ p 10 N 1000 P 10 12 15 20 25 30 35 40"

#tmux new-session -d  "Rscript aggregate_results.R path simulationsNEW/ p 20 P 20"
#tmux new-session -d  "Rscript aggregate_results.R path simulationsNEW/ p 30 P 30"
#tmux new-session -d  "Rscript aggregate_results.R path simulationsNEW/ p 40 P 40"
#tmux new-session -d  "Rscript aggregate_results.R path simulationsNEW/ p 50 P 50"
