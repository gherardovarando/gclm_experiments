#!/bin/bash
tmux new-session -d  "Rscript simulate_ou.R path simulations/ p 10 P 10"
tmux new-session -d  "Rscript simulate_ou.R path simulations/ p 10 P 12"
tmux new-session -d  "Rscript simulate_ou.R path simulations/ p 10 P 15"
tmux new-session -d  "Rscript simulate_ou.R path simulations/ p 10 P 20"
tmux new-session -d  "Rscript simulate_ou.R path simulations/ p 10 P 25"
tmux new-session -d  "Rscript simulate_ou.R path simulations/ p 10 P 30"
tmux new-session -d  "Rscript simulate_ou.R path simulations/ p 10 P 35"
tmux new-session -d  "Rscript simulate_ou.R path simulations/ p 10 P 40"

tmux new-session -d  "Rscript simulate_ou.R path simulations/ p 20 P 20"
tmux new-session -d  "Rscript simulate_ou.R path simulations/ p 30 P 30"
tmux new-session -d  "Rscript simulate_ou.R path simulations/ p 40 P 40"
tmux new-session -d  "Rscript simulate_ou.R path simulations/ p 50 P 50"
