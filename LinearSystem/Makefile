CFLAGS = -O2 -Wall -Wfatal-errors

.PHONY: all


all: LinearSystem.h5 LinearSystem.so LinearSystem_gvf.c

LinearSystem_gvf.c: LinearSystem.vf
	vfgen gsl LinearSystem.vf

LinearSystem.so: LinearSystem_gvf.c LinearSystem_gvf.h
	gcc -shared $(CFLAGS) -o $@ $<

LinearSystem.h5: CreateData.sh LinearSystem.vf InitialValue.txt time.txt trajectory.double
	./CreateData.sh

LinearSystem.vf: MakeLinearSystemUsingGNUOctave.m
	./MakeLinearSystemUsingGNUOctave.m

