#!/bin/bash
gcc -Wall -c bsm.c -std=c99
g++ -Wall -c Prob_bsm.cc
g++ -Wall Prob_bsm.o bsm.o -lglobes -lgsl -lgslcblas -O3 -o Prob_bsm
./Prob_bsm 
rm *.o
