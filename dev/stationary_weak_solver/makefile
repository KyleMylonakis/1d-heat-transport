OBJS = functions.cpp dynamics.cpp Noise.cpp Weights_and_Nodes.cpp randomize_initial_data.cpp
CC = g++
CFLAGS = -Wall -std=gnu++17 -fopenmp -lfftw3 -lm -O3

dynamics : 
	$(CC) $(OBJS) $(CFLAGS) -o dynamics

clean: 
	\rm *.o *~ dynamics
