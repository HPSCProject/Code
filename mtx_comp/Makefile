CC = icc
LD = icc

OPT = -Ofast
mtx_cmp.exe: mtx_cmp.o 
	$(LD) $(OPT) -o mtx_cmp.exe mtx_cmp.o

mtx_cmp.o: mtx_cmp.cpp 
	$(CC) $(OPT) -c mtx_cmp.cpp 


clean:
	rm -f *.o Breathing.exe 
