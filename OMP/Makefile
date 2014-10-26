CC = icc
FC = ifort
LD = ifort

OPT = -Ofast -openmp -vec_report6 -simd -xhost

Breathing.exe: Breathing.o functions.o 
 	$(CC) $(OPT) -o Breathing.exe Breathing.o functions.o 

Breathing.o: Breathing.cpp 
	$(CC) $(OPT) -c Breathing.cpp 


functions.o: functions.cpp 
	$(CC) $(OPT) -c functions.cpp 

clean:
	rm -f *.o Breathing.exe 
