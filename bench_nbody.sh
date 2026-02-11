#!/bin/bash
#SBATCH --job-name=nbody_bench
#SBATCH --error=nbody_%j.err
#SBATCH --time=05:00:00
#SBATCH --partition=Centaurus
#SBATCH --mem=30G

srun $HOME/Nobody/nbody solar.tsv 200 5000000 10000 solar_out1.tsv
srun $HOME/Nobody/nbody 100 1 10000 100 part_out1.tsv
srun $HOME/Nobody/nbody 1000 1 10000 100 part_out2.tsv
