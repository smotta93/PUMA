#!/bin/bash
gcc -Wall -c bsm.c -std=c99
g++ -Wall -c Prob_bsm_loop.cc
g++ -Wall Prob_bsm_loop.o bsm.o -lglobes -lgsl -lgslcblas -O3 -o Prob_bsm_loop
./Prob_bsm_loop 
rm *.o
