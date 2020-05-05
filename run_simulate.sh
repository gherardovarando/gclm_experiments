#!/bin/bash
tmux new-session -d  "Rscript simulate.R path simulations/ p 10 N 1000 P 10 "
tmux new-session -d  "Rscript simulate.R path simulations/ p 10 N 1000 P 12 "
tmux new-session -d  "Rscript simulate.R path simulations/ p 10 N 1000 P 15 "
tmux new-session -d  "Rscript simulate.R path simulations/ p 10 N 1000 P 20 "
tmux new-session -d  "Rscript simulate.R path simulations/ p 10 N 1000 P 25 "
tmux new-session -d  "Rscript simulate.R path simulations/ p 10 N 1000 P 30 "
tmux new-session -d  "Rscript simulate.R path simulations/ p 20 N 1000 P 20 "
tmux new-session -d  "Rscript simulate.R path simulations/ p 30 N 1000 P 30 "
tmux new-session -d  "Rscript simulate.R path simulations/ p 40 N 1000 P 40 "
tmux new-session -d  "Rscript simulate.R path simulations/ p 50 N 1000 P 50 "
tmux new-session -d  "Rscript simulate.R path simulations/ p 100 N 1000 P 100 "

