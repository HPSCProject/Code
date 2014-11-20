#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=12
#SBATCH --time=00:30:00
#SBATCH --job-name=jasmine_project
#SBATCH --output=project.out

./8-16_breathing.exe 1001 1001 50 500 >> out.txt

