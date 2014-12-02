CC = icc
LD = icc

OPT = -Ofast -openmp -simd -xhost -vec_report6

all: 11-25_serial.exe 11-25_breathing_omp.exe

11-25_serial.exe: 11-25_serial.o
	$(LD) $(OPT) -o 11-25_serial.exe 11-25_serial.o

11-25_breathing.exe: 11-25_breathing.o
	$(LD) $(OPT) -o 11-25_breathing.exe 11-25_breathing.o

11-25_breathing_omp.exe: 11-25_breathing_omp.o
	$(LD) $(OPT) -o 11-25_breathing_omp.exe 11-25_breathing_omp.o

#Breathing.exe: Breathing.o functions.o 
#	$(LD) $(OPT) -o Breathing.exe Breathing.o functions.o 

#Breathing.o: Breathing.cpp 
#	$(CC) $(OPT) -c Breathing.cpp 

#functions.o: functions.cpp 
#	$(CC) $(OPT) -c functions.cpp

11-25_serial.o:11-25_serial.cpp
	$(CC) $(OPT) -c 11-25_serial.cpp

11-25_breathing.o:11-25_breathing.cpp
	$(CC) $(OPT) -c 11-25_breathing.cpp

11-25_breathing_omp.o:11-25_breathing_omp.cpp
	$(CC) $(OPT) -c 11-25_breathing_omp.cpp

run:
	./11-25_breathing.exe 125 125 125 25.0 250 

clean:
	rm -f *.o 11-25_breathing.exe 11-25_breathing_omp.exe 11-25_serial.exe 
