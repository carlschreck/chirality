#!/bin/bash 

# Remove opengl program if "make" doesn't recognize changes
rm opengl

# Compile imaging program
make

# Run for example
./opengl ../program/prod_v9_ar1.3_L20.0_layer13.0_desync0.4_seed2.0_b4e3_twist15.0.dat 20.0
