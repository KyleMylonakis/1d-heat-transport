OBJS = functions.cpp dynamics.cpp Noise.cpp
CC = g++
CFLAGS = -Wall -std=gnu++1z -fopenmp -lfftw3 -lm -O3

dynamics : 
	$(CC) $(OBJS) $(CFLAGS) -o dynamics

clean: 
	\rm *.o *~ dynamics
