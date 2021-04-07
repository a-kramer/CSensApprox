CC = gcc
CFLAGS = --std=gnu11 -march=native -O2 -Wall -Wfatal-errors `pkg-config --cflags gsl blas hdf5`
LDFLAGS = `pkg-config --libs gsl blas hdf5` -lhdf5_hl 
.PHONY: all


all: gsl_odeiv


h5block.o: h5block.c h5block.h
	gcc $(CFLAGS) -c $< $(LDFLAGS)

gsl_odeiv: main.c h5block.o ndarray.o solution.o
	gcc $(CFLAGS) -o $@ $^ $(LDFLAGS) -ldl


