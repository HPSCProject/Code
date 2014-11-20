CC = icc
LD = icc

OPT = -Ofast -vec_report6 -simd -xhost

Breathing.exe: Breathing.o
	$(LD) -o Breathing.exe Breathing.o 

Breathing.o: Breathing.cpp 
	$(CC) $(OPT) -c Breathing.cpp 

clean:
	rm -f *.o Breathing.exe 
