#!/bin/bash

## Specifies the interpreting shell for the job.
#$ -S /bin/bash

## Specifies that all environment variables active within the qsub utility be exported to the context of the job.
#$ -V

## Specifies the parallel environment
#$ -pe smp 4

## Execute the job from the current working directory.
#$ -cwd 

## The  name  of  the  job.
#$ -N MandelbrotOpenMP

##send an email when the job ends
#$ -m e

##email addrees notification
#$ -M rtv3@alumnes.udl.cat

##Passes an environment variable to the job
#$ -v  OMP_NUM_THREADS=16

## In this line you have to write the command that will execute your application.
./mandelbrot_openmp 600 400 10000 
