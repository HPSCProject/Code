CC = icc
LD = icc


OPT = -Ofast
Breathing.exe: Breathing.o functions.o 
	$(LD) $(OPT) -o Breathing.exe Breathing.o functions.o

Breathing.o: Breathing.cpp 
	$(CC) $(OPT) -c Breathing.cpp 

functions.o: functions.cpp 
	$(CC) $(OPT) -c functions.cpp 

clean:
	rm -f *.o Breathing.exe 
