CC = icc

OPT = -Ofast

24-10_D3Q41.exe: timing.o 24-10_D3Q41.o
	$(CC) $(OPT) -o 24-10_D3Q41.exe timing.o 24-10_D3Q41.o

24-10_D3Q41.o: 24-10_D3Q41.cpp
	$(CC) $(OPT) -c 24-10_D3Q41.cpp

timing.o: timing.c
	$(CC) $(OPT) -c timing.c
