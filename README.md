# Parallel-Mandelbrot
This project aims to solve a specific problem by means of MPI and OpenMP
programming models. Specifically, the main aim is to understand the given real
problem, to provide different parallel solutions and analyse the scalability of the
provided solutions and the effects of the design decisions.

## Description of the problem
We want to implement a parallel version for a popular program, the Mandelbrot set. The
Mandelbrot set is a geometric figure of infinite complexity (fractal nature) obtained from
a mathematical formula and a small recursive algorithm, as shown in the following
figure:

<p align="center">
  <img src="https://github.com/Rochii/OpenMP-Mandelbrot/blob/master/serial/mandelbrot_serial.png?raw=true" width="650" alt="Mandelbrot Set Output">
</p>

The code provided in this project is a serial implementation (in C) of the Mandelbrot set
that is shown in the previous Figure, where the intensity of the background colour
indicates the proximity of each point to an element of the set (Black: element belongs to
the set, White: it is close to the set and Blue: far from the set). 

## OpenMP implementation
In order to obtain the best performance of an application implemented with a
parallel programming model, we should start by optimizing the application at
node-level. According to this, the first activity will be focused on study the
operation of the sequential version, analysing the pieces of code candidate to be
parallelized following a work and data decomposition model. Next, we will apply
the OpenMP directives to parallelize the corresponding code according to the
hardware and the suitable work decomposition model at node level.

