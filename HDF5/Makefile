CC = icc

INCLUDE = -I$(H5DIR)/include
LIB_DIR = -L$(H5DIR)/lib
LIBS = -lhdf5

8-16_breathing:
	$(CC) 8-16_breathing.cpp $(INCLUDE) $(LIB_DIR) $(LIBS) -o 8-16_breathing.exe

clean:
	rm *.o *.exe
