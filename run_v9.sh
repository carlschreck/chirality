#! /bin/bash -l

#SBATCH --qos=shared
#SBATCH --constraint=haswell
#SBATCH --time=48:00:00
#SBATCH --ntasks=1

ar=$1
Lx=$2
D1=$3

twist=$4

layerwidth=$5
layerdepth=$6
frontdepth=$7

propdepth=$8
bounddepth=$9
traildepth=${10}

rate0=${11}
desync=${12}
seed=-${13}

numsteps=${14}
dataskip=${15}
prodskip=${16}
layerskip=${17}
dt=${18}
b=${19}

movie=${20}

prodfile=prod_v9_ar${ar}_L${Lx}_layer${layerdepth}_desync${desync}_seed$7_b${b}_twist${twist}

rundir=~/production/linear_front_ellipse_damped_removepart_twist
outdir=/global/cscratch1/sd/cschreck/linear_front_ellipse_damped_removepart_twist

cd $outdir

time $rundir/ellipse_divide_linear_damped_twist_v9.o << EOF
  $ar
  $Lx
  $D1
  $twist
  $layerwidth
  $layerdepth
  $frontdepth
  $propdepth
  $bounddepth
  $traildepth
  $rate0
  $desync
  $seed
  $numsteps
  $dataskip
  $prodskip
  $layerskip
  $dt
  $b
  $movie
  $prodfile
EOF
