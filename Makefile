CC = gcc
CFLAGS = --std=gnu11 -march=native -O2 -Wall -Wfatal-errors
.PHONY: all


all: gsl_odeiv


h5block.o: h5block.c h5block.h
	gcc $(CFLAGS) -c `pkg-config --cflags gsl hdf5`  $^ `pkg-config --libs gsl hdf5` -lhdf5_hl 


gsl_odeiv: main.c h5block.o
	gcc $(CFLAGS) `pkg-config --cflags gsl hdf5` -o $@ $^ `pkg-config --libs gsl hdf5` -lhdf5_hl -ldl


