#!/bin/bash


if [ -z $2 ] ;
then
	echo "please specify two worldfiles!"
else
	echo "Starting Simulation with worldfile: $1"
	
	#echo "Running: ./Polyworld $1"
	./Polyworld $1
	
	mkdir -p jay_runs/
	mkdir -p jay_runs/$1
	
	if [ -z jay_runs/$1/run ];
	then
		mv run jay_runs/$1/
	else
		rm -rf jay_runs/$1/run
		mv run jay_runs/$1/
	fi
	

	echo "Generating seedGenomes.txt and seedPositions.txt"
	./scripts/genomeSeed --pos jay_runs/$1/run
	
	ls ./genomeSeeds.txt
	ls ./seedPositions.txt
	
	echo "Starting transplant simulation using worldfile: $2"
	
	./Polyworld $2
	
	
fi



