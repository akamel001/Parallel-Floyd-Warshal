#!/bin/bash

#script that submits jobs and variees the threads
#while keeping number of nodes constant 

for i in {1..12}
do
	qsub run-omp.qsub $i 500
done

for i in {1..12}
do
	qsub run-omp.qsub $i 1000
done

for i in {1..12}
do
	qsub run-omp.qsub $i 2000
done

for i in {1..12}
do
	qsub run-omp.qsub $i 3000
done
