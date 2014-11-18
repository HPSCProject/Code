#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=12
#SBATCH --time=05:00:00
#SBATCH --job-name=jasmine_project
#SBATCH --output=project.out

./24-10_D3Q41.exe 501 501 501 100 250 >> out.txt
