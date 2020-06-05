#include <stdio.h>
#include <complex.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "mpi.h"
#include <omp.h>

typedef unsigned char pixel_t[3]; // colors [R, G ,B]

// Main program
int main(int argc, char *argv[])
{
    int x, y;                                           /* Each iteration, it calculates: newz = oldz*oldz + p, where p is the current pixel, and oldz stars at the origin */
    int rank, namelen, size, num;                       /* MPI world variables */
    int W, H, MAXITER;                                  /* Image width, height and number of mandel iterations */
    double pr, pi;                                      /* Real and imaginary part of the pixel p */    
    double newRe, newIm, oldRe, oldIm;                  /* Real and imaginary parts of new and old z */    
    double zoom = 1, moveX = -0.5, moveY = 0;           /* You can change these to zoom and change position */
    double tic, toc, time_s, time_e;                    /* MPI time variables */
    char host[50];                                      /* Host buffer */
    int *myRecvArr = (int *)malloc(H*W*sizeof(pixel_t) + 1);   /* Receive MPI array */
    int *mySendArr = (int *)malloc(H*W*sizeof(pixel_t) + 1);   /* Send MPI array */

    MPI_Status status;                                  /* MPI status variable */

    // Argument parsing
    if(argc == 4)
    {
        W = atoi(argv[1]);
        H = atoi(argv[2]);
        MAXITER = atoi(argv[3]);
    }
    else{
        printf(" Invalid parameters. Usage: mandelbrot_openmp <WIDTH> <HEIGHT> <ITERATIONS>\n");
        exit(0);
    }

    printf("Executing mandelbrot with widht: %d, height: %d for %d iterations.\n", W, H, MAXITER);

    // MPI initilizations
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);           /* Get the rank of the process */
    MPI_Comm_size(MPI_COMM_WORLD, &size);           /* Get all the processes */
    MPI_Get_processor_name(host, &namelen);         /* Get the host */
    time_s = MPI_Wtime();                           
    tic = clock();


    if(rank == 0) // Master
    {
        // Get the number of chunks (number of rows)
        int n_chunks = 0;
        int total_chunks = H;

        mySendArr[0] = n_chunks;
        // While we have chunks
        while(n_chunks < total_chunks)
        {
            // Receive worker requests
            MPI_Recv(myRecvArr, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
            printf(" MASTER => received request from worker[%d]\n", status.MPI_SOURCE);

            // Worker want more work
            if(myRecvArr[0] == -1)
            {
                // Send chunk to the worker and increase the next chunk to send
                MPI_Send(mySendArr, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
                printf(" MASTER => received work request from worker[%d]\n", status.MPI_SOURCE);

                n_chunks++;
                mySendArr[0] = n_chunks;
            }
            // Worker want to send finished work
            else if(myRecvArr[0] == -2)
            {
                // Receive pixels from workers
                MPI_Recv(myRecvArr, H*W + 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD, &status);
                printf(" MASTER => received request to get worker[%d] results\n", status.MPI_SOURCE);
                
                // TODO: Append data to somewhere

                // TODO: Clean receive array
                free(myRecvArr);
            }
        }

        time_e = MPI_Wtime();
        MPI_Finalize();
        toc = clock();
        printf("Task: %d - Time (Wtime): %f - Time (clock): %f\n", rank, (time_e-time_s), (double)((toc-tic) / CLOCKS_PER_SEC));        
    }
    else // Worker
    {
        int i;

        while(1)
        {
            // Ask master for work
            mySendArr[0] = -1;
            MPI_Send(mySendArr, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
            printf(" Process[%d] => sending work request to master\n", rank);

            // Recv response (starting number or -1)
            MPI_Recv(myRecvArr, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
            printf(" Process[%d] => received master message\n", rank);

            // -1 means no more data to process, break the loop
            if(myRecvArr[0] == -1){
                printf(" Process[%d] => breaking form the loop\n", rank);
                break;
            }
            else
            {
                pixel_t *pixels = malloc(sizeof(pixel_t)*H*W + 1);  /* Reserve memory to allocate colour pixels */

                y = myRecvArr[0];
                printf(" Process[%d] => completing chunk %d\n", rank, y);

                #pragma omp parallel for shared(pixels, moveX, moveY, zoom) private(x, y, pr, pi, newRe, newIm, oldRe, oldIm)  schedule(dynamic)
                /* loop through every pixel */                
                for(x = 0; x < W; x++)
                {
                    /* Calculate the initial real and imaginary part of z, based on the pixel location and zoom and position values */
                    pr = 1.5 * (x - W/2) / (0.5*zoom*W) + moveX;
                    pi = (y - H/2) / (0.5*zoom*H) + moveY;
                    newRe = newIm = oldRe = oldIm = 0;
                    
                    /* Start the iteration process */
                    for(i = 0; i < MAXITER; i++)
                    {
                        /* Remember value of previous iteration */
                        oldRe = newRe;
                        oldIm = newIm;

                        /* The actual iteration, the real and imaginary part are calculated */
                        newRe = oldRe * oldRe - oldIm * oldIm + pr;
                        newIm = 2 * oldRe * oldIm + pi;
                        
                        /* If the point is outside the circle with radius 2: stop */
                        if((newRe * newRe + newIm * newIm) > 4) break;
                    }
                    
                    if(i == MAXITER)
                    {
                        pixels[y*W + x][0] = 0;
                        pixels[y*W + x][1] = 0;
                        pixels[y*W + x][2] = 0;
                    }
                    else
                    {
                        double z = sqrt(newRe*newRe + newIm*newIm);
                        int brightness = 256 * log2(1.75 + i - log2(log2(z))) / log2((double)MAXITER);
                        pixels[y*W + x][0] = brightness;
                        pixels[y*W + x][1] = brightness;
                        pixels[y*W + x][2] = 255;                
                    }    
                }

                printf(" Process[%d] => finished chunk %d calculation\n", rank, y);

                // Tell master work is done and ready to send
                mySendArr[0] = -2;
                MPI_Send(mySendArr, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
                printf(" Process[%d] => sending that I finished the work with chunk %d\n", rank, y);

                // Send work
                mySendArr[0] = y; // Chunk row
                /*for(i = 1; i < y*W+1; i++){
                    mySendArr[i] = pixels[i-1];
                }*/

                // Send all pixels
                //MPI_Send(mySendArr, 1*W + 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
                printf(" Process[%d] => Sending all pixels calculed\n", rank);

                //free(pixels);
            }            
        }
        MPI_Finalize();
    }

    return 0;
}

/*
        // Write pixels in a file
        char filename[50];
        int y_act, x_act;

        sprintf(filename, "../files/mandelbrot_hybrid_%d.ppm", rank);
        FILE *fp = fopen(filename, "wb");
        fprintf(fp, "P6\n# CREATOR: Roger Truchero\n");
        fprintf(fp, "%d %d\n255\n", W, (end-begin));
        
        for(y_act = begin; y_act < end; y_act++){                
            for(x_act = 0; x_act < W; x_act++){
                fwrite(pixels[y_act*W + x_act], 1, sizeof(pixel_t), fp);
            }
        }

        fclose(fp);
*/