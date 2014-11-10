#!/bin/bash
#SBATCH -N 1                      # Number of requested nodes
#SBATCH --ntasks-per-node=12      # number of cores per node
#SBATCH -t 1:30:00                # Requests a maximum CPU time of three minutes
#SBATCH -J job               # Names the job
#SBATCH --output=job.out     # Output file name
#SBATCH --qos janus         # Select the debug queue

./11-2_breathing.exe 251 251 50 250
./11-2_breathing.exe 501 501 100 250
./11-2_breathing.exe 1001 1001 200 250

./serial_breathing.exe 251 251 50 250
./serial_breathing.exe 501 501 100 250
./serial_breathing.exe 1001 1001 200 250
 

