CFLAGS = -O2 -Wall -Wfatal-errors

.PHONY: all

all: CaMKIIs.so

%.so: %_gvf.c %_gvf.h
	gcc -shared $(CFLAGS) -o $@ $<

