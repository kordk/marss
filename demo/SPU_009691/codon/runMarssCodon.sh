#!/bin/sh -x

TWOBIT=../strPur4.local.1kb.2bit
GPRED=../SPU_009691.local.ref.gp
BAMLIST=spp_to_bam.region.txt

time -p ~/proj/marss/src/marssCodon \
    -n 8 \
    -2 $TWOBIT \
    -p $GPRED \
    -b $BAMLIST \
    -m freq -D -g > marssCodon.out 2>&1

