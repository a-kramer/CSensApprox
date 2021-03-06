#!/bin/bash

H5F=LinearSystem.h5
DATASET=/data/trajectory

if [[ -f LinearSystem.h5 ]]; then
    rm LinearSystem.h5
fi

h5import trajectory.double -d 6,$((2**7)) -p $DATASET -s 64 -o $H5F
h5import mu.double -d 6 -p /prior/mu -s 64 -o $H5F

if (($?==0)); then
    ../h5attr -d $DATASET -a index -s 0 $H5F
    ../h5attr -d $DATASET -a time time.txt $H5F
    ../h5attr -d $DATASET -a input -s 1 $H5F
    ../h5attr -d $DATASET -a InitialValue InitialValue.txt $H5F
fi
