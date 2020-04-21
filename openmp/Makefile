all:
	gcc mandelbrot_openmp.c -o mandelbrot_openmp -fopenmp -lm
	./mandelbrot_openmp
	convert mandelbrot_openmp.ppm mandelbrot_openmp.png 
	xdg-open mandelbrot_openmp.png

clean:
	rm mandelbrot_openmp
	rm mandelbrot_openmp.png
