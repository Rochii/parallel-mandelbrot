all:
	gcc mandelbrot_openmp.c -o mandelbrot_openmp -fopenmp -lm
	./mandelbrot_openmp
	convert mandelbrot_openmp.ppm mandelbrot_openmp.png
	rm mandelbrot_openmp.ppm
	xdg-open mandelbrot_openmp.png
