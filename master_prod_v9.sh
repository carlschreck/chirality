#!/bin/bash

ar=3.0
Lx=5.0
D1=1.0

twist=15.0

layerwidth=1.0
layerdepth=13.0
frontdepth=2.0

propdepth=4.0
bounddepth=6.0
traildepth=2.0

rate0=1.0
desync=0.4

numsteps=25000000 # longer run for production
dataskip=100000
prodskip=100000
layerskip=200
dt=4e-6
b=4e3

movie=.TRUE.

for seed in `seq 101 110`;
do
  sbatch run_v9.sh $ar $Lx $D1 $twist $layerwidth $layerdepth $frontdepth $propdepth $bounddepth $traildepth $rate0 $desync $seed $numsteps $dataskip $prodskip $layerskip $dt $b $movie
done
