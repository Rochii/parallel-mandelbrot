static:
	mpicc mandelbrot_hybrid_static.c -o mandelbrot_hybrid_static -fopenmp -lm

dynamic: 
	mpicc mandelbrot_hybrid_dynamic.c -o mandelbrot_hybrid_dynamic -fopenmp -lm

run: 
	qsub run_script.sh

clean: 
	rm output error ../files/*
