#!/usr/bin/bash

mpicxx main.cpp --std=c++11 -o mTaichi
mpirun -np 2 mTaichi
rm mTaichi
exit

