# Twist program

Fortran program for linear front of ellipse-shaped cells that twist while growing

## Important parameters

There are a number of parameters in this simulation. The most important are:

twist - twist rate (degrees / cell cycle)

layerdepth - growth layer depth

alpha0 - aspect ratio at birth

b - damping coefficient# chirality

The codes contained in this repo are:

1. ellipse_divide_linear_damped_twist_v9.f - Self-contained Fortran program that models cellular growth via proliferation of ellipse-shaped cells and overdamped molecular dynamics. 

2. ellipse_divide_linear_damped_twist_v9.o - Executable for ellipse_divide_linear_damped_twist_v9.f

3. run_v9.sh - based script to run ellipse_divide_linear_damped_twist_v9.o

4. master_example_v9.sh - Example bash script to run run_v9.sh. The purpose of this script is to loop over parameters of interest. 

5. master_prod_v9.sh - Production bash scriot to submit jobs to slurm queue. This script loops over intial random seed.
