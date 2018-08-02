OBJS = functions.cpp dynamics.cpp
CC = g++
CFLAGS = -Wall -std=gnu++1z -fopenmp -O3

dynamics : 
	$(CC) $(OBJS) $(CFLAGS) -o dynamics

clean: 
	\rm *.o *~ dynamics
