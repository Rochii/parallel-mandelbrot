all:
	mpicc mandelbrot_mp.c -o mandelbrot_mp -lm

run: 
	qsub run_script.sh

clean: 
	rm output error ../files/*
