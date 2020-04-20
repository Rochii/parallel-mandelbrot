
/* Includes*/
# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <omp.h>

/* Functions */
int i4_min(int i1, int i2);

/* Defines */
#define M 500
#define N 500
#define ITER 2000

/* Main program */
int main()
{    
    int r[M][N];
    int g[M][N];
    int b[M][N];   
    int count[M][N];  
    int c, i, j, jhi, jlo, k;

    char *output_filename = "mandelbrot_openmp.ppm";
    FILE *output_unit;

    double wtime;

    double x_max =   1.25;
    double x_min = - 2.25;
    double x;
    double x1;
    double x2;
    double y_max =   1.75;
    double y_min = - 1.75;
    double y;
    double y1;
    double y2;


    printf("Starting mandelbrot execution!\n");

    wtime = omp_get_wtime ( );

    /* Carry out the iteration for each pixel, determining COUNT */
    # pragma omp parallel shared (b, count, g, r, x_max, x_min, y_max, y_min) private (i, j, k, x, x1, x2, y, y1, y2)
    {
        # pragma omp for
        for(i = 0; i < M; i++)
        {
            y = ((double)(i-1)*y_max + (double)(M-i)*y_min) / (double)(M-1);

            for(j = 0; j < N; j++)
            {
                x = ((double)(j-1)*x_max + (double)(N-j)*x_min) / (double)(N-1);
                count[i][j] = 0;

                x1 = x;
                y1 = y;

                for(k = 1; k <= ITER; k++)
                {
                    x2 = x1*x1 - y1*y1 + x;
                    y2 = 2*x1*y1 + y;

                    if(x2 < -2.0 || 2.0 < x2 || y2 < -2.0 || 2.0 < y2)
                    {
                        count[i][j] = k;
                        break;
                    }
                    x1 = x2;
                    y1 = y2;
                }

                if(k >= ITER){
                    r[i][j] = 0; g[i][j] = 0; b[i][j] = 0;
                }
                else
                {
                    c = (int)(255.0 * sqrt(sqrt(sqrt(((double)(count[i][j]) / (double)(ITER))))));
                    r[i][j] = 3 * c / 5; 
                    g[i][j] = 3 * c / 5; 
                    b[i][j] = c;
                }
            }
        }
    }

    wtime = omp_get_wtime() - wtime;
    printf("\n  Time = %g seconds.\n", wtime);
    
    /* Write data to an ASCII PPM file. */
    output_unit = fopen(output_filename, "wt");
    fprintf(output_unit, "P3\n");
    fprintf(output_unit, "%d  %d\n", N, M);
    fprintf(output_unit, "%d\n", 255);

    for(i = 0; i < M; i++)
    {
        for(jlo = 0; jlo < N; jlo = jlo + 4)
        {
            jhi = jlo+4 < N ? jlo+4 : N;

            for(j = jlo; j < jhi; j++){
                fprintf(output_unit, "  %d  %d  %d", r[i][j], g[i][j], b[i][j]);
            }
            fprintf(output_unit, "\n");
        }
    }
    /* Close file descriptor */
    fclose(output_unit);
    printf("\n  Graphics data written to \"%s\".\n", output_filename);

    return 0;
}