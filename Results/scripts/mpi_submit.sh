#!/bin/bash

#mpi script fixed nodes at 750 
NODES=400

for i in {1..12} 
do
	qsub run-mpi.qsub $i 750
done


for j in {1..12} 
do
	qsub run-mpi2.qsub 7 $NODES
	let NODES=NODES+200
done

