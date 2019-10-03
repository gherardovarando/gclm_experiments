#!/bin/bash
tmux new-session -d -s 11 "Rscript simulate_ou.R path simulations2/ p 10 P 10"
tmux new-session -d -s 12 "Rscript simulate_ou.R path simulations2/ p 10 P 12"
tmux new-session -d -s 13 "Rscript simulate_ou.R path simulations2/ p 10 P 15"
tmux new-session -d -s 14 "Rscript simulate_ou.R path simulations2/ p 10 P 20"
tmux new-session -d -s 15 "Rscript simulate_ou.R path simulations2/ p 10 P 25"
tmux new-session -d -s 16 "Rscript simulate_ou.R path simulations2/ p 10 P 30"
tmux new-session -d -s 17 "Rscript simulate_ou.R path simulations2/ p 10 P 35"
tmux new-session -d -s 18 "Rscript simulate_ou.R path simulations2/ p 10 P 40"
