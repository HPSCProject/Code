CC = icc

OPT = -Ofast

8-16_breathing.exe: timing.o 8-16_breathing.o
	$(CC) $(OPT) -o 8-16_breathing.exe timing.o 8-16_breathing.o

8-16_breathing.o: 8-16_breathing.cpp
	$(CC) $(OPT) -c 8-16_breathing.cpp

timing.o: timing.c
	$(CC) $(OPT) -c timing.c
