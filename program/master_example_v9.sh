#!/bin/bash

gfortran -O3 -o ellipse_divide_linear_damped_twist_v9.o ellipse_divide_linear_damped_twist_v9.f

ar=2.0
D1=1.6
Lx=100.0

twist=100.0

layerdepth=13.0
layerwidth=1.0
frontdepth=2.0

propdepth=4.0
bounddepth=6.0
traildepth=2.0

rate0=0.2
desync=0.4
seed=101

numsteps=35000000
dataskip=100000
prodskip=20000 
layerskip=100
dt=4e-6
b=4e3
bani=0.0

movie=.TRUE.

./run_v9.sh $ar $Lx $D1 $twist $layerwidth $layerdepth $frontdepth $propdepth $bounddepth $traildepth $rate0 $desync $seed $numsteps $dataskip $prodskip $layerskip $dt $b $bani $movie
