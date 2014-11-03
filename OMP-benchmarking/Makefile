CC = icc

OPT = -Ofast -openmp -simd #-vec_report6

all: serial_breathing.exe 11-2_breathing.exe #8-16_breathing.exe

11-2_breathing.exe: 11-2_breathing.o
	$(CC) $(OPT) -o 11-2_breathing.exe 11-2_breathing.o

#8-16_breathing.exe: timing.o 8-16_breathing.o
#	$(CC) $(OPT) -o 8-16_breathing.exe timing.o 8-16_breathing.o

serial_breathing.exe: serial_breathing.o
	$(CC) $(OPT) -o serial_breathing.exe serial_breathing.o  

11-2_breathing.o: 11-2_breathing.cpp
	$(CC) $(OPT) -c 11-2_breathing.cpp

#8-16_breathing.o: 8-16_breathing.cpp
#	$(CC) $(OPT) -c 8-16_breathing.cpp

serial_breathing.o: serial_breathing.cpp
	$(CC) $(OPT) -c serial_breathing.cpp

#timing.o: timing.c
#	$(CC) $(OPT) -c timing.c

run:
#	./11-2_breathing.exe 251 251 50 250
	./11-2_breathing.exe 501 501 50 250
#	./11-2_breathing.exe 1001 1001 200 250
#	./serial_breathing.exe 251 251 50 250
	./serial_breathing.exe 501 501 50 250
#	./serial_breathing.exe 1001 1001 200 250
clean:
	rm -f *.o serial_breathing.exe 11-2_breathing.exe #8-16_breathing.exe
