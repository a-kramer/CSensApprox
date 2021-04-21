CFLAGS = -O2 -Wall -Wfatal-errors

.PHONY: all

all: HarmonicOscillator.so

%.so: %_gvf.c %_gvf.h
	gcc -shared $(CFLAGS) -o $@ $<

HarmonicOscillator.h5: MakeHarmonicOscillator.R
	./MakeHarmonicOscillator.R

