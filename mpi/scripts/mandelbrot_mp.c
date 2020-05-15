#include <stdio.h>
#include <complex.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "mpi.h"

#define SIZE 800

// Function definitions
int mandel(complex z0);

// Main program
int main(int argc, char *argv[])
{
    double xmin, xmax, ymin, ymax;
    int i, j, rows, columns;
    complex z;
    int row[SIZE];
    unsigned char line[3 * SIZE];
    FILE *img;
    char file[80];

    // MPI variables initilizations
    int rank, namelen, size, num;
    long acc = 0;
    double result, tic, toc, time_s, time_e;
    char host[50], cont;
    MPI_Status status;  

    // MPI initilizations
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Get_processor_name(host,&namelen);
    time_s = MPI_Wtime();
    tic = clock();

    int n = (SIZE * SIZE) / (size - 1);
    xmin = -2; xmax = 1;
    ymin = -1.5; ymax = 1.5;
  
    if(rank == 0) // Master
    {
        snprintf(file, sizeof file, "%s%s", "mandel", ".pam");
        img = fopen(file,"w");
        fprintf(img, "P6\n%d %d 255\n", SIZE, SIZE);
        for(i = 1; i <= size - 1; i++) 
        {
            int begin = ((i - 1) * n) / SIZE;
            int end = ((i * n - 1)) / SIZE;
            for(j = begin; j <= end; j++)
            {
                MPI_Recv(&line, 3 * SIZE, MPI_UNSIGNED_CHAR, i, 1, MPI_COMM_WORLD, &status);
                fwrite(line, 1, sizeof line, img);	
            }
        }
    }
    else
    {   // Worker 	
        int begin = ((rank - 1) * n) / SIZE;
        int end = ((rank * n - 1)) / SIZE;

        for(i = begin; i <= end; i++)
        {
            for(j = 0; j < SIZE; j++)
            {
                z = xmin + j*((xmax-xmin)/SIZE) + (ymax - i*((ymax-ymin)/SIZE))*I;
                row[j] = mandel(z);
            }

            for(j = 0; j < SIZE; j++)
            {
                if(row[j] <= 63)
                {
                    line[3*j] = 255;
                    line[3*j+1] = line[3*j + 2] = 255 - 4*row[j];
                }
                else
                {
                    line[3*j] = 255;
                    line[3*j + 1] = row[j] - 63;
                    line[3*j + 2] = 0;
                }

                if(row[j] == 320) line[3*j] = line[3*j + 1] = line[3*j + 2] = 255;
            }
            MPI_Send(&line, 3 * SIZE, MPI_UNSIGNED_CHAR, 0, 1, MPI_COMM_WORLD);
        }
    }

    time_e = MPI_Wtime();
    MPI_Finalize();
    toc = clock();
    printf("Task: %d - Time (Wtime): %ld - Time (clock): %ld\n", rank, (time_e-time_s), (double)((toc-tic) / CLOCKS_PER_SEC));

    return 0;
}


int mandel(complex z0)
{
    int i;
    complex z;

    z = z0;
    for(i = 1; i < 320; i++)
    {
        z = z*z + z0;
        if((creal(z)*creal(z)) + (cimag(z)*cimag(z)) > 4.0){
            break;
        }
    }
    
    return i;
}