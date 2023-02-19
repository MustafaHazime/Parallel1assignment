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

    if (rank == 0) {
        // Master process
        int num_workers = size - 1;
        int next_row = 0;

        while (next_row < HEIGHT) {
            // Send the next task to a worker
            MPI_Status status;
            int worker_rank;
            MPI_Recv(&worker_rank, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            int row = next_row;
            next_row++;

            if (row >= HEIGHT) {
                // No more tasks to do, send termination signal to worker
                int tag = 0;
                MPI_Send(&row, 1, MPI_INT, worker_rank, tag, MPI_COMM_WORLD);
            } else {
                // Send the task to the worker
                int tag = 1;
                MPI_Send(&row, 1, MPI_INT, worker_rank, tag, MPI_COMM_WORLD);
            }
        }
    } else {
        // Worker process
        while (1) {
            // Receive the next task from the master
            int row;
            MPI_Status status;
            int tag = MPI_ANY_TAG;
            MPI_Recv(&row, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);

            if (tag == 0) {
                // Termination signal, exit the loop
                break;
            } else {
                // Compute the Mandelbrot set for this row
                int image_row[WIDTH];
                for (int x = 0; x < WIDTH; x++) {
                    double real = X_MIN + (double)x / (WIDTH - 1) * (X_MAX - X_MIN);
                    double imag = Y_MIN + (double)row / (HEIGHT - 1) * (Y_MAX - Y_MIN);
                    int iterations = mandelbrot(real, imag, MAX_ITERATIONS);
                    image_row[x] = iterations;
                }

                // Send the result back to the master
                MPI_Send(image_row, WIDTH, MPI_INT, 0, row, MPI_COMM_WORLD);
            }
        }
    }

    MPI_Finalize();

end = clock();
time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
printf("Time taken: %f seconds\n", time_used);
    return 0;
}
