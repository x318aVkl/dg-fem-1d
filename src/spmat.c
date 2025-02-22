

#include "spmat.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>



#define SPMAT_UNSET (INT_MAX-1)


sparse_matrix_t sparse_matrix_create(
    int n_rows,
    int n_per_row
) {
    sparse_matrix_t mat;


    mat.n_rows = n_rows;
    mat.n_per_row = n_per_row;

    mat.columns = (int*)malloc(sizeof(int) * n_rows * n_per_row);
    mat.values = (double*)malloc(sizeof(double) * n_rows * n_per_row);

    // set all columns to unallocated and values to zero
    for (int i=0; i<(n_rows * n_per_row); ++i) {
        mat.columns[i] = SPMAT_UNSET;
        mat.values[i] = 0.0;
    }

    return mat;
}


void sparse_matrix_free(
    sparse_matrix_t spmat
) {
    free(spmat.columns);
    free(spmat.values);
}



double* sparse_matrix_index(
    sparse_matrix_t spmat,
    int i,
    int j
) {
    int start = i * spmat.n_per_row;
    int end = (i + 1) * spmat.n_per_row;

    if ((i < 0) || (i >= spmat.n_rows)) {
        fprintf(stderr, "ERROR in sparse_matrix_index, index is %d > matrix size %d\n", i, spmat.n_rows);
        exit(1);
    }

    int k = start;
    while ((spmat.columns[k] != j) && (spmat.columns[k] != SPMAT_UNSET) && (k < end)) {
        k++;
    }

    if (k == end) {
        fprintf(stderr, "ERROR in sparse_matrix_index, row %d is full\n", i);
        exit(1);
    }

    if (spmat.columns[k] == SPMAT_UNSET) {
        spmat.columns[k] = j;
    }

    return spmat.values + k;
}



int gauss_siedel(
    sparse_matrix_t lhs,
    const double* rhs,
    double* x,
    double omega,
    double tol,
    int verbose
) {

    int n = lhs.n_rows;

    // init x to zero
    for (int i=0; i<n; ++i)
        x[i] = 0.0;

    double residual = 1.0;

    double norm_rhs = 0.0;
    for (int i=0; i<n; ++i) {
        norm_rhs += rhs[i] * rhs[i];
    }
    norm_rhs = sqrt(norm_rhs);

    int iter = 0;
    while (residual > tol) {

        residual = 0.0;

        for (int it=0; it<n; ++it) {
            int i;
            if (iter % 2 == 0) {
                i = it;
            }  else {
                i = (n - 1) - it;
            }

            double sigma = 0.0;
            double a_ii = 0.0;
            for (int k = i*lhs.n_per_row; k < (i+1)*lhs.n_per_row; ++k) {
                int j = lhs.columns[k];
                if (j == SPMAT_UNSET) {
                    break;
                }

                if (i == j) {
                    a_ii = lhs.values[k];
                } else {
                    sigma += lhs.values[k] * x[j];
                }
            }
            
            if (a_ii < 1e-14) {
                fprintf(stderr, "ERROR in gauss_siedel, zero diagonal\n");
                exit(1);
            }

            double r_i = sigma + a_ii * x[i] - rhs[i];

            double x_i = (rhs[i] - sigma) / a_ii;
            double dx = x_i - x[i];

            x[i] = x[i] * (1.0 - omega) + x_i * omega;


            if (verbose > 1) {
                printf("    %d %lf\n", i, x[i]);
            }
            
            residual += r_i * r_i;
        }

        residual = sqrt(residual) / norm_rhs;

        if (verbose > 0) {
            printf("%d %.3e\n", iter, residual);
        }

        iter++;
    }

    return iter;
}



