CFLAGS = -O2 -Wall -Wfatal-errors

.PHONY: all clean


all: LinearSystem.h5 LinearSystem.so LinearSystem_gvf.c


clean:
	rm LinearSystem.h5 LinearSystem.so LinearSystem_gvf.[ch]

LinearSystem_gvf.c: LinearSystem.vf
	vfgen gsl LinearSystem.vf

LinearSystem.so: LinearSystem_gvf.c LinearSystem_gvf.h
	gcc -shared -fPIC $(CFLAGS) -o $@ $<

LinearSystem.h5: CreateData.sh LinearSystem.vf InitialValue.txt time.txt trajectory.double mu.double
	./CreateData.sh

#LinearSystem.vf: MakeLinearSystemUsingGNUOctave.m
#	./MakeLinearSystemUsingGNUOctave.m

