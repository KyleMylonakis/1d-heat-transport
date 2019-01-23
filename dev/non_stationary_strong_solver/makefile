OBJS = functions.cpp dynamics.cpp Weights_and_Nodes.cpp
CC = g++
CFLAGS = -Wall -std=gnu++17 -fopenmp -O3

dynamics : 
	$(CC) $(OBJS) $(CFLAGS) -o dynamics

clean: 
	\rm *.o *~ dynamics
