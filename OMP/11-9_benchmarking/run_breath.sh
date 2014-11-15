#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=12
#SBATCH --time=00:30:00
#SBATCH --job-name=jasmine_hw05
#SBATCH --output=breath.out

export OMP_NUM_THREADS="2"
./Breathing.exe >> 1001.txt

export OMP_NUM_THREADS="4"
./Breathing.exe >> 1001.txt

export OMP_NUM_THREADS="6"
./Breathing.exe >> 1001.txt

export OMP_NUM_THREADS="8"
./Breathing.exe >> 1001.txt

export OMP_NUM_THREADS="10"
./Breathing.exe >> 1001.txt

export OMP_NUM_THREADS="12"
./Breathing.exe >> 1001.txt
