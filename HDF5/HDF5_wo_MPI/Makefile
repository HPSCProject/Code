CC = icc

INCLUDE = -I$(H5DIR)/include
LIB_DIR = -L$(H5DIR)/lib
LIBS = -lhdf5

HDF5_Example:
	$(CC) HDF5_Example.cpp $(INCLUDE) $(LIB_DIR) $(LIBS) -o HDF5_Example.exe

clean:
	rm *.o *.exe
