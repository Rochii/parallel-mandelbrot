all:
	mpicc mandelbrot_hybrid.c -o mandelbrot_hybrid -fopenmp -lm

run: 
	qsub run_script.sh

clean: 
	rm output error ../files/*
