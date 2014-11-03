#!/bin/bash

#SBATCH -N 1		# Requests one node
#SBATCH -n 12		# Requests 12 cores per node
#SBATCH -t 3:00:00	# Requests max CPU time of 3 hrs
#SBATCH -J Benchmarking	# Assigns name to job
#SBATCH --qos normal	# Assign to a normal Janus node


export OMP_NUM_THREADS=2
echo "Number of threads = 2" > benchmarking.out
./11-2_breathing.exe 251 251 50 250 >> benchmarking.out
./11-2_breathing.exe 251 251 50 500 >> benchmarking.out
./11-2_breathing.exe 501 501 50 250 >> benchmarking.out
./11-2_breathing.exe 501 501 50 500 >> benchmarking.out
./11-2_breathing.exe 1001 1001 50 250 >> benchmarking.out
./11-2_breathing.exe 1001 1001 50 500 >> benchmarking.out

export OMP_NUM_THREADS=4
echo "Number of threads = 4" >> benchmarking.out
./11-2_breathing.exe 251 251 50 250 >> benchmarking.out
./11-2_breathing.exe 251 251 50 500 >> benchmarking.out
./11-2_breathing.exe 501 501 50 250 >> benchmarking.out
./11-2_breathing.exe 501 501 50 500 >> benchmarking.out
./11-2_breathing.exe 1001 1001 50 250 >> benchmarking.out
./11-2_breathing.exe 1001 1001 50 500 >> benchmarking.out

export OMP_NUM_THREADS=6
echo "Number of threads = 6" >> benchmarking.out
./11-2_breathing.exe 251 251 50 250 >> benchmarking.out
./11-2_breathing.exe 251 251 50 500 >> benchmarking.out
./11-2_breathing.exe 501 501 50 250 >> benchmarking.out
./11-2_breathing.exe 501 501 50 500 >> benchmarking.out
./11-2_breathing.exe 1001 1001 50 250 >> benchmarking.out
./11-2_breathing.exe 1001 1001 50 500 >> benchmarking.out

export OMP_NUM_THREADS=8
echo "Number of threads = 8" >> benchmarking.out
./11-2_breathing.exe 251 251 50 250 >> benchmarking.out
./11-2_breathing.exe 251 251 50 500 >> benchmarking.out
./11-2_breathing.exe 501 501 50 250 >> benchmarking.out
./11-2_breathing.exe 501 501 50 500 >> benchmarking.out
./11-2_breathing.exe 1001 1001 50 250 >> benchmarking.out
./11-2_breathing.exe 1001 1001 50 500 >> benchmarking.out

export OMP_NUM_THREADS=10
echo "Number of threads = 10" >> benchmarking.out
./11-2_breathing.exe 251 251 50 250 >> benchmarking.out
./11-2_breathing.exe 251 251 50 500 >> benchmarking.out
./11-2_breathing.exe 501 501 50 250 >> benchmarking.out
./11-2_breathing.exe 501 501 50 500 >> benchmarking.out
./11-2_breathing.exe 1001 1001 50 250 >> benchmarking.out
./11-2_breathing.exe 1001 1001 50 500 >> benchmarking.out

export OMP_NUM_THREADS=12
echo "Number of threads = 12" >> benchmarking.out
./11-2_breathing.exe 251 251 50 250 >> benchmarking.out
./11-2_breathing.exe 251 251 50 500 >> benchmarking.out
./11-2_breathing.exe 501 501 50 250 >> benchmarking.out
./11-2_breathing.exe 501 501 50 500 >> benchmarking.out
./11-2_breathing.exe 1001 1001 50 250 >> benchmarking.out
./11-2_breathing.exe 1001 1001 50 500 >> benchmarking.out

echo "OMP Benchmarking is finished" > message.txt
/bin/mail -s "OMP BENCHMARKING" "samuel.elliott@colorado.edu" < message.txt
rm message.txt
