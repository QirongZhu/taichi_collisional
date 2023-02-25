#!/usr/bin/bash

mpicxx main.cpp -Wall --std=c++11 -o mTaichi -O3 -mtune=native
mpirun -np 2 mTaichi
rm mTaichi
exit

