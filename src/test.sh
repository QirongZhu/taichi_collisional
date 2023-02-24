#!/usr/bin/bash

mpicxx main.cpp -Wall --std=c++11 -o mTaichi 
mpirun -np 2 mTaichi
rm mTaichi
exit

