#!/bin/bash
#SBATCH --time=3
#SBATCH --nodes=1
#SBATCH --mem=1GB
#$ -cwd
#$ -S /bin/bash


reps=1
ns=(25 75)
B=50000
theta=.2

## students t call
for n in ${ns[@]}; do for x in $(seq 1 $reps); do for param in $(seq 2.1 .1 5); do Rscript sim.R $B $n 1 $param $theta; done; done& done


## power law call:
for n in ${ns[@]}; do for x in $(seq 1 $reps); do for param in $(seq -.9 .1 0); do Rscript sim.R $B $n 2 $param $theta; done; done &done

## beta call
for n in ${ns[@]}; do for x in $(seq 1 $reps); do for param in $(seq .1 .1 1); do Rscript sim.R $B $n 3 $param $theta; done; done &done
