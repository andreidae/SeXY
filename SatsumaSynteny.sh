#!/bin/bash

#1 - Threads
#2 - Query sequence (genome of interest)
#3 - Target sequence (Chromosomes)
#4 - Output_directory

/groups/hologenomics/andreaac/Software/satsuma/SatsumaSynteny -m 1 -n $1 -q $2 -t $3 -o $4
