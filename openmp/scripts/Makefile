all: 
	gcc mandelbrot_openmp.c -o mandelbrot_openmp -fopenmp -lm

run:
	./mandelbrot_openmp 600 400 10000
	convert ../files/mandelbrot_openmp.ppm ../files/mandelbrot_openmp.png 
	xdg-open ../files/mandelbrot_openmp.png

clean:
	rm ../files/mandelbrot_openmp.ppm
	rm ../files/mandelbrot_openmp.png
