CC=gcc

# CFLAGS will be the options passed to the compiler. 
CFLAGS= -O3 -std=c++11 
#CFLAGS= -g -std=c++11 -fopenmp 
OBJECTS  = main.o functions.o 

all: 2D_HOCFD

2D_HOCFD:  $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o $@ 

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< 

clean: 
	rm -rf   *.o 
