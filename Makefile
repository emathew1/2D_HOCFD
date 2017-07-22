CC=g++

# CFLAGS will be the options passed to the compiler. 
CFLAGS= -O3 -std=c++11 
#CFLAGS= -g -std=c++11 -fopenmp 
OBJECTS  = main.o Utils.o Solver.o Filter.o Derivatives.o IdealGas.o 

all: 2D_HOCFD

2D_HOCFD:  $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o $@ 

main.o: main.cpp Utils.hpp Filter.hpp Derivatives.hpp IdealGas.hpp
	$(CC) $(CFLAGS) -c $< 

Utils.o: Utils.cpp Utils.hpp
	$(CC) $(CFLAGS) -c $<

Solver.o: Solver.cpp Solver.hpp Utils.hpp Filter.hpp Derivatives.hpp IdealGas.hpp
	$(CC) $(CFLAGS) -c $<

Filter.o: Filter.cpp Filter.hpp Utils.hpp
	$(CC) $(CFLAGS) -c $<

Derivatives.o: Derivatives.cpp Derivatives.hpp Utils.hpp
	$(CC) $(CFLAGS) -c $<

IdealGas.o: IdealGas.cpp 
	$(CC) $(CFLAGS) -c $<

clean: 
	rm -rf   *.o 
