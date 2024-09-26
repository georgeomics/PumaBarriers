#!/bin/bash
#author: ge199066
#SBATCH -J tess3r
#SBATCH -o tess3r.out
#SBATCH -p normal
#SBATCH --cpus-per-task=16
#SBATCH -t 7-00:00:00
#SBATCH --mem-per-cpu=16000

cd /home/ge199066/projects/puma-AZ/pop-structure/TESS/coombs/SCRIPTS

module load R/431

Rscript tess3r_coombs.R

