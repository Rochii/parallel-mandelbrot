#include <stdio.h>
#include <complex.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "mpi.h"
#include <omp.h>

//typedef unsigned char pixel_t[3]; // colors [R, G ,B]

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
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);               /* Get the rank of the process */
    MPI_Comm_size(MPI_COMM_WORLD, &size);               /* Get all the processes */
    MPI_Get_processor_name(host, &namelen);             /* Get the host */
    time_s = MPI_Wtime();                           
    tic = clock();

    int myRecvArr[W+1][3];                                 /* Receive MPI array */
    int mySendArr[W+1][3];                                 /* Send MPI array */

    // Initialize matrix's with 0's
    for(x = 0; x < W+1; x++)
    {
        for(y = 0; y < 3; y++){
            myRecvArr[x][y] = 0;
            mySendArr[x][y] = 0;            
        }
    }

    if(rank == 0) // Master
    {
        // Get the number of chunks (number of rows)
        int n_chunks = 0;
        int total_chunks = H;
        int recv_chunk = 0;
        int recv_chunks = 0;
        int workers_end = 0;
        int pixel_result[H][W][3];        
        
        mySendArr[0][0] = n_chunks;

        // Until we have received all the pieces, still waiting for requests
        while(recv_chunks < total_chunks || workers_end < size-1)
        {
            printf(" MASTER => waiting for worker request. n_chunks: %d total_chunks: %d recv_chunks: %d workers_end: %d\n", n_chunks, total_chunks, recv_chunks, workers_end);        

            // Receive worker requests
            MPI_Recv(&myRecvArr, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
            printf(" MASTER => received request from worker[%d]\n", status.MPI_SOURCE);

            // Worker want more work
            if(myRecvArr[0][0] == -1)
            {
                // If no more work to send
                if(n_chunks == total_chunks)
                {
                    printf(" MASTER => no more work to send, n_chunks: %d\n", n_chunks);
                    mySendArr[0][0] = -1;
                    MPI_Send(&mySendArr, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
                    workers_end++;
                }
                else
                {              
                    mySendArr[0][0] = n_chunks;

                    // Send chunk to the worker and increase the next chunk to send
                    printf(" MASTER => sending chunk %d to worker[%d]\n", n_chunks, status.MPI_SOURCE);
                    MPI_Send(&mySendArr, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);

                    n_chunks++;                    
                    printf(" MASTER => Next chunk to send: %d\n", n_chunks);
                }                
            }
            // Worker want to send finished work
            else if(myRecvArr[0][0] == -2)
            {
                // Receive pixels from workers
                MPI_Recv(&myRecvArr, W+1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD, &status);
                printf(" MASTER => received request to get worker[%d] results\n", status.MPI_SOURCE);

                recv_chunk = myRecvArr[0][0];
                recv_chunks++;

                printf(" MASTER => Received chunk: %d Received chunks: %d\n", recv_chunk, recv_chunks);              
                                
/*
                // TODO: Append data to pixel_result
                for(x = 1; x < W+1; x++)
                {
                    for(y = 0; y < 3; y++){
                        pixel_result[recv_chunk][x][y];
                    }
                }
*/
            }
        }

        printf(" MASTER => All chunks received\n");

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
            mySendArr[0][0] = -1;
            MPI_Send(&mySendArr, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
            printf(" Process[%d] => sending work request to master\n", rank);

            // Recv response (starting number or -1)
            MPI_Recv(&myRecvArr, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
            printf(" Process[%d] => received master message\n", rank);

            // -1 means no more data to process, break the loop
            if(myRecvArr[0][0] == -1){
                printf(" Process[%d] => breaking form the loop\n", rank);
                break;
            }
            else
            {
                int pixels[W][3];  /* Reserve memory to allocate colour pixels */

                y = myRecvArr[0][0];
                printf(" Process[%d] => completing chunk %d\n", rank, y);

                #pragma omp parallel for shared(pixels, moveX, moveY, zoom) private(x, y, pr, pi, newRe, newIm, oldRe, oldIm)  schedule(dynamic)
                /* loop through every pixel */                
                for(x = 0; x < W; x++)
                {
                    //printf(" Process[%d] => x: %d\n", rank, x);

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
                        newRe = oldRe*oldRe - oldIm*oldIm + pr;
                        newIm = 2*oldRe*oldIm + pi;
                        
                        /* If the point is outside the circle with radius 2: stop */
                        if((newRe * newRe + newIm * newIm) > 4) break;
                    }
                    
                    if(i == MAXITER)
                    {
                        pixels[x][0] = 0;
                        pixels[x][1] = 0;
                        pixels[x][2] = 0;
                    }
                    else
                    {
                        double z = sqrt(newRe*newRe + newIm*newIm);
                        int brightness = 256 * log2(1.75 + i - log2(log2(z))) / log2((double)MAXITER);
                        pixels[x][0] = brightness;
                        pixels[x][1] = brightness;
                        pixels[x][2] = 255;  
                    }
                }

                printf(" Process[%d] => finished chunk %d calculation\n", rank, y);

                // Tell master work is done and ready to send
                mySendArr[0][0] = -2;
                MPI_Send(&mySendArr, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
                printf(" Process[%d] => sending that I finished the work with chunk %d\n", rank, y);

                // Send work
                mySendArr[0][0] = y; // Chunk row
                int k;
                for(i = 1; i < W+1; i++){
                    for(k = 0; k < 3; k++){
                        mySendArr[i][k] = pixels[i-1][k];
                    }
                }

                // Send all pixels to the MASTER
                MPI_Send(&mySendArr, W+1, MPI_INT, 0, 0, MPI_COMM_WORLD);
                printf(" Process[%d] => Sending all pixels calculated\n", rank);
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