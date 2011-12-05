#!/bin/bash

#script that submits jobs and variees the threads
#while keeping number of nodes constant 
NODES=400

for i in {1..12}
do
	qsub run-omp.qsub 			$i 3000
	qsub run-fwomp.qsub 			$i 3000
	qsub run-mpi.qsub 			$i 3000 
	qsub run-mpi-complex.qsub	$i 3000
done

for j in {1..12} 
do
	qsub run-omp.qsub 			8 $NODES
	qsub run-fwomp.qsub 			8 $NODES
	qsub run-mpi.qsub 			8 $NODES
	qsub run-mpi-complex.qsub	8 $NODES
	let NODES=NODES+300
done
