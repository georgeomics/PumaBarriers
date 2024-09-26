#!/bin/bash
#author: ge199066
#SBATCH -J structure
#SBATCH -o structure.out
#SBATCH -e structure.err
#SBATCH -p normal
#SBATCH --cpus-per-task=16
#SBATCH -t 7-00:00:00
#SBATCH --mem-per-cpu=16000

cd /home/ge199066/projects/puma-AZ/pop-structure/STRUCTURE/coombs/DATA

module load structure

for i in {1..10};
do
  for j in {1..10};
  do
    structure -K "$i" -o outputfile_K"$i"_"$j"
  done
done
