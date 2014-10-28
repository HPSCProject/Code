CC = icc

OPT = -Ofast -openmp

8-16_breathing.exe: timing.o 8-16_breathing.o
	$(CC) $(OPT) -o 8-16_breathing.exe timing.o 8-16_breathing.o

8-16_breathing.o: 8-16_breathing.cpp
	$(CC) $(OPT) -c 8-16_breathing.cpp

timing.o: timing.c
	$(CC) $(OPT) -c timing.c

run:
	./8-16_breathing.exe 501 501 100 250

clean:
	rm -f *.o 8-16_breathing.exe 
