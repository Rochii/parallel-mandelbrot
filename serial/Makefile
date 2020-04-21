all:
	gcc mandelbrot_serial.c -o mandelbrot_serial -lm
	./mandelbrot_serial > mandelbrot_serial.ppm
	convert mandelbrot_serial.ppm mandelbrot_serial.png
	xdg-open mandelbrot_serial.png

clean:
	rm mandelbrot_serial
	rm mandelbrot_serial.png
