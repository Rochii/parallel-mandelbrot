all:
	gcc mandelbrot_openmp.c -o mandelbrot_openmp -fopenmp -lm
	./mandelbrot_openmp
	convert files/mandelbrot_openmp.ppm files/mandelbrot_openmp.png 
	xdg-open files/mandelbrot_openmp.png

clean:
	rm files/mandelbrot_openmp
	rm files/mandelbrot_openmp.png
