all:
	mpicc mandelbrot_mpi.c -o mandelbrot_mpi -lm

run: 
	qsub run_script.sh

clean: 
	rm output error ../files/*
