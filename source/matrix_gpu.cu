/* ****************************************************************************************************************** */
/**
 *  @file matrix_gpu.cu
 *  @author Edward J. Parkinson (e.parkinson@soton.ac.uk)
 *  @date August 2023
 *
 *  @brief
 *
 *  ***************************************************************************************************************** */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#if CUDA_ON

#include <cuda.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>

/* NVCC is a C++ compiler at heart, so anything we re-use from regular C source
   has to be defined here with `extern "C"` to tell the compiler that the
   function has been compiled by a C compiler (and does some computer science
   stuff to make linking possible) */

extern "C" int Exit (int error_code);
extern "C" int Error (const char *format, ...);
extern "C" int Log (const char *format, ...);

/* `cusolver_handle` is a variable used to interact with the cuSolver/CUDA
    runtime and is used to initialise and clean up the resources required for
    both runtimes */

static cusolverDnHandle_t cusolver_handle = NULL;

/* ****************************************************************************************************************** */
/**
 *  @brief Check the return status of a CUDA function
 *
 *  @param [in] status  the status to check
 *
 *  ***************************************************************************************************************** */

#define CUDA_CHECK(status)                                                                                             \
  do {                                                                                                                 \
    cudaError_t err = status;                                                                                          \
    if (err != cudaSuccess) {                                                                                          \
      Error("CUDA Error: %s\n", cudaGetErrorString(err));                                                              \
      Exit(EXIT_FAILURE);                                                                                              \
    }                                                                                                                  \
  } while (0)

/* ****************************************************************************************************************** */
/**
 *  @brief Check the return status of a cuSOLVER function
 *
 *  @param [in] status  the status to check
 *
 *  ***************************************************************************************************************** */

#define CUSOLVER_CHECK(status)                                                                                         \
  do {                                                                                                                 \
    cusolverStatus_t err = status;                                                                                     \
    if (err != CUSOLVER_STATUS_SUCCESS) {                                                                              \
      Error("cuSolver Error: %d\n", err);                                                                              \
      Exit(EXIT_FAILURE);                                                                                              \
    }                                                                                                                  \
  } while (0)

/* ****************************************************************************************************************** */
/**
 * @brief
 *
 * @return
 *
 * @details
 *
 *  ***************************************************************************************************************** */

extern "C" void
create_cusolver_handle (void)
{
  CUSOLVER_CHECK (cusolverDnCreate (&cusolver_handle));
  Log ("Created a new cuSOLVER handle created\n");
}

/* ****************************************************************************************************************** */
/**
 * @brief
 *
 * @return
 *
 * @details
 *
 *  ***************************************************************************************************************** */

extern "C" void
destroy_cusolver_handle (void)
{
  CUSOLVER_CHECK (cusolverDnDestroy (cusolver_handle));
  Log ("Destroyed the cuSOLVER handle\n");
}

/* ****************************************************************************************************************** */
/**
 * @brief  Solve the linear system A x = b, for the vector x
 *
 * @param  [in]  a_matrix  a square matrix on the LHS
 * @param  [in]  b_vector  the B resultant vector
 * @param  [in]  size  the number of rows (and columns) in the square matrix matrix and vectors
 * @param  [out] x_vector  the x vector on the RHS
 *
 * @return an integer representing the error state
 *
 * @details
 * Performs LU decomposition to solve for x in the linear system A x = b. The calculation is perform, in serial, on
 * the CPU using GSL.
 *
 *  ***************************************************************************************************************** */

extern "C" int
gpu_solve_linear_system (double *a_matrix, double *b_vector, int size, double *x_vector)
{
  if (cusolver_handle == NULL)
    create_cusolver_handle ();

  int *devInfo;
  int lwork;
  int *d_pivot;                 /* device array of pivoting sequence */
  double *d_A, *d_b;
  double *d_work;               /* cuSolver needs a "workspace" to do stuff, which we have to allocate manually */

  /* Allocate memory on the GPU (device) to store the matrices/vectors */
  cudaMalloc ((void **) &d_A, size * size * sizeof (double));
  cudaMalloc ((void **) &d_b, size * sizeof (double));
  cudaMalloc ((void **) &devInfo, sizeof (int));

  /* Copy the matrix and vector to the device memory */
  cudaMemcpy (d_A, a_matrix, size * size * sizeof (double), cudaMemcpyHostToDevice);
  cudaMemcpy (d_b, b_vector, size * sizeof (double), cudaMemcpyHostToDevice);
  cudaMalloc ((void **) &d_pivot, size * sizeof (int));

  /* XXXX_bufferSize is used to compute the size of the workspace we need, and depends on the size of the linear
     system being solved */
  cusolverDnDgetrf_bufferSize (cusolver_handle, size, size, d_A, size, &lwork);
  cudaMalloc ((void **) &d_work, lwork * sizeof (double));

  /* Perform LU factorization and solve the linear system. The vector d_b is not used in `getrs` (the solver), but
     it's the same size of the solution vector so we'll re-use that. d_b is then copied back to host memory (CPU RAM) */
  cusolverDnDgetrf (cusolver_handle, size, size, d_A, size, d_work, d_pivot, devInfo);
  cusolverDnDgetrs (cusolver_handle, CUBLAS_OP_T, size, 1, d_A, size, d_pivot, d_b, size, devInfo);
  cudaMemcpy (x_vector, d_b, size * sizeof (double), cudaMemcpyDeviceToHost);

  cudaFree (d_A);
  cudaFree (d_b);
  cudaFree (d_work);
  cudaFree (d_pivot);

  return EXIT_SUCCESS;
}

/* ****************************************************************************************************************** */
/**
 * @brief
 *
 * @param  [in]  a_matrix
 * @param  [out] a_inverse
 * @param  [in]  num_rows
 *
 * @return an integer representing the error state
 *
 * @details
 *
 *  ***************************************************************************************************************** */

extern "C" int
gpu_invert_matrix (double *matrix, double *inverse, int num_rows)
{
  if (!cusolver_handle)
    create_cusolver_handle ();

  return EXIT_SUCCESS;
}

#endif
