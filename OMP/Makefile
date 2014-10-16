CC = icc
FC = ifort
LD = ifort

OPT = -Ofast -openmp -vec_report6 -simd -xhost

Breathing.exe: Breathing.o functions.o 
 	$(LD) -o Breathing.exe Breathing.o functions.o 

Breathing.o: Breathing.cpp 
	$(FC) $(OPT) -c Breathing.cpp 


functions.o: functions.cpp 
	$(FC) $(OPT) -c functions.cpp 

clean:
	rm -f *.o Breathing.exe 
