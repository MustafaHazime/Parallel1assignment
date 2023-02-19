#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>

int mandelbrot(double x, double y, int max_iterations) {
    double real = 0.0;
    double imag = 0.0;
    int iterations = 0;
    while (real * real + imag * imag < 4.0 && iterations < max_iterations) {
        double temp_real = real * real - imag * imag + x;
        double temp_imag = 2.0 * real * imag + y;
        real = temp_real;
        imag = temp_imag;
        iterations++;
    }
    return iterations;
}

int main(int argc, char **argv) {
  clock_t start, end;
  double time_used;
  start = clock();
    const int WIDTH = 800;
    const int HEIGHT = 600;
    const double X_MIN = -2.0;
    const double X_MAX = 1.0;
    const double Y_MIN = -1.5;
    const double Y_MAX = 1.5;
    const int MAX_ITERATIONS = 1000;

    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int rows_per_process = HEIGHT / size;
    int start_row = rank * rows_per_process;
    int end_row = start_row + rows_per_process - 1;

    // Handle remainder rows
    if (rank == size - 1) {
        end_row = HEIGHT - 1;
    }

    int image_rows[end_row - start_row + 1][WIDTH];

    // Compute the Mandelbrot set for this process's rows
    for (int row = start_row; row <= end_row; row++) {
        for (int x = 0; x < WIDTH; x++) {
            double real = X_MIN + (double)x / (WIDTH - 1) * (X_MAX - X_MIN);
            double imag = Y_MIN + (double)row / (HEIGHT - 1) * (Y_MAX - Y_MIN);
            int iterations = mandelbrot(real, imag, MAX_ITERATIONS);
            image_rows[row - start_row][x] = iterations;
        }
    }

    // Gather the results from all processes
    if (rank == 0) {
        int result[HEIGHT][WIDTH];
        MPI_Gather(image_rows, (end_row - start_row + 1) * WIDTH, MPI_INT, result, (end_row - start_row + 1) * WIDTH, MPI_INT, 0, MPI_COMM_WORLD);

        // Print the Mandelbrot set
        for (int row = 0; row < HEIGHT; row++) {
            for (int x = 0; x < WIDTH; x++) {
                if (result[row][x] == MAX_ITERATIONS) {
                    printf("*");
                } else {
                    printf(" ");
                }
            }
            printf("\n");
        }
    } else {
        MPI_Gather(image_rows, (end_row - start_row + 1) * WIDTH, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
    }

    MPI_Finalize();
end = clock();
time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
printf("Time taken: %f seconds\n", time_used);

    return 0;
}
