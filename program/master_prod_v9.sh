#!/bin/bash
  
ar=1.3
D1=1.6
Lx=20.0

twist=15.0

layerdepth=13.0
layerwidth=1.0
frontdepth=2.0

propdepth=4.0
bounddepth=6.0
traildepth=2.0

rate0=1.0
desync=0.4
seed=101 

numsteps=25000000
dataskip=100000
prodskip=100000
layerskip=100
dt=4e-6
b=4e3
    
movie=.TRUE.

for seed in `seq 1 10`;
do
  sbatch run_v9.sh $ar $Lx $D1 $twist $layerwidth $layerdepth $frontdepth $propdepth $bounddepth $traildepth $rate0 $desync $seed $numsteps $dataskip $prodskip $layerskip $dt $b $movie
done

