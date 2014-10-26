CC = mpicxx
LD = mpicxx

OPT = -Ofast -vec_report6 -simd -xhost

Breathing.exe: Breathing_2d.o functions.o 
	$(LD) -o Breathing.exe Breathing_2d.o functions.o 

Breathing_2d.o: Breathing_2d.cpp 
	$(CC) $(OPT) -c Breathing_2d.cpp 

functions.o: functions.cpp 
	$(CC) $(OPT) -c functions.cpp 

clean:
	rm -f *.o Breathing.exe 
