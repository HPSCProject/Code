CC = mpicxx
LD = mpicxx

OPT = -O3 -openmp

Breathing.exe: Breathing_2d.o functions.o 
	$(LD) $(OPT) -o Breathing.exe Breathing_2d.o functions.o 

Breathing_2d.o: Breathing_2d.cpp 
	$(CC) $(OPT) -c Breathing_2d.cpp -o Breathing_2d.o

functions.o: functions.cpp 
	$(CC) $(OPT) -c functions.cpp -o functions.o

clean:
	rm -f *.o Breathing.exe 
