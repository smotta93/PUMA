#!/bin/bash

# Compilar o código
g++ -Wall -c Prob_std.cc
g++ -Wall Prob_std.o -lglobes -lgsl -lgslcblas -O3 -o Prob_std

# Executar o programa
./Prob_std

# Remover arquivos objeto
rm *.o
