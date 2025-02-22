
#include "dg_assemble.h"
#include "spmat.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


// problem selection
//#define PROBLEM_POISSON
#define PROBLEM_CONVECT


#define MAX_DEGREE 12


// forward declaration of assembly functions
double interior_lhs(shape_t u, shape_t v, double dx, assembly_data_t data, void* user_data);
double interior_rhs(shape_t v, double dx, assembly_data_t data, void* user_data);
double numerical_flux(
    shape_t u0, shape_t v0,
    shape_t u1, shape_t v1,
    double n, double ds, double h, 
    assembly_data_t data,
    void* user_data
);
double boundary_term_lhs(
    shape_t u, shape_t v,
    double n, double ds, double h, 
    assembly_data_t data,
    void* user_data
);
double boundary_term_rhs(
    shape_t v,
    double n, double ds, double h, 
    assembly_data_t data,
    void* user_data
);



// main code
int main(int argc, char** argv) {

    if (argc < 3) {
        printf("example <size> <degree>\n");
        exit(1);
    }

    int size = atoi(argv[1]);
    int degree = atoi(argv[2]);

    if ((size < 1) || (size > 100000)) {
        fprintf(stderr, "error, invalid size %d\n", size);
        exit(1);
    }
    if ((degree < 0) || (degree > 100)) {
        fprintf(stderr, "error, invalid degree %d\n", degree);
        exit(1);
    }

    printf("Running with size %d, degree %d\n", size, degree);


    double dx = 1.0 / ((double)size);

    int n_dofs = size * (degree + 1);
    int dofs_per_cell = degree + 1;

    printf("- n_dofs = %d\n", n_dofs);
    
    sparse_matrix_t full_lhs = sparse_matrix_create(
        n_dofs,
        dofs_per_cell * 3
    );
    
    double* full_rhs = (double*)malloc(sizeof(double) * n_dofs);
    for (int i=0; i<n_dofs; ++i)
        full_rhs[i] = 0.0;    

    // pack boundary functions into arrays
    assembly_boundary_lhs_t boundary_lhs[] = {
        boundary_term_lhs
    };
    assembly_boundary_rhs_t boundary_rhs[] = {
        boundary_term_rhs
    };

    // make the basis
    //basis_t basis = legendre_basis(degree);
    basis_t basis = lagrange_basis(degree);

    printf("- Assembling linear system\n");
    // loop over cells and assemble problem
    for (int cell = 0; cell < size; ++cell) {

        // create local systems
        double local_lhs_diag[MAX_DEGREE*MAX_DEGREE];
        double local_lhs_left[MAX_DEGREE*MAX_DEGREE];
        double local_lhs_right[MAX_DEGREE*MAX_DEGREE];
        double local_rhs[MAX_DEGREE];

        double* local_lhs[] = {
            local_lhs_diag,
            local_lhs_left,
            local_lhs_right
        };

        // create boundary information
        int bound_left = -1;
        if (cell == 0) bound_left = 0;

        int bound_right = -1;
        if (cell == (size - 1)) bound_right = 0;
        
        int boundaries[] = {bound_left, bound_right};

        // assemble the local system
        assemble_local_system(
            local_lhs,
            local_rhs,
            basis,
            dx,
            boundaries,
            cell,
            interior_lhs,
            interior_rhs,
            numerical_flux,
            boundary_lhs,
            boundary_rhs,
            (void*)&degree
        );

        // distribute the local system
        for (int i=0; i<dofs_per_cell; ++i) {
            full_rhs[cell * dofs_per_cell + i] += local_rhs[i];
        }

        for (int i=0; i<dofs_per_cell; ++i) {
            int gi = cell * dofs_per_cell + i;
            for (int j=0; j<dofs_per_cell; ++j) {

                int gj0 = cell * dofs_per_cell + j;

                *sparse_matrix_index(full_lhs, gi, gj0) += local_lhs[0][i*dofs_per_cell + j];
                
                if (cell > 0) {
                    int gj1 = (cell - 1) * dofs_per_cell + j;

                    *sparse_matrix_index(full_lhs, gi, gj1) += local_lhs[1][i*dofs_per_cell + j];
                }
                if (cell < (size - 1)) {
                    int gj2 = (cell + 1) * dofs_per_cell + j;

                    *sparse_matrix_index(full_lhs, gi, gj2) += local_lhs[2][i*dofs_per_cell + j];
                }
            }
        }
    }



    // solve the problem
    double* u = (double*)malloc(sizeof(double) * n_dofs);
    int iterations = gauss_siedel(
        full_lhs, full_rhs, u, 
        // change omega convergence parameter for the matrix
        #if defined(PROBLEM_POISSON)
            1.0, 
        #elif defined(PROBLEM_CONVECT)
            0.1,
        #else
            #error "Define PROBLEM_POISSON or PROBLEM_CONVECT"
        #endif
        1e-12, 0
    );
    printf("- Solved problem in %d iterations\n", iterations);

    // compute the l2 norm of the error
    double l2_error = 0.0;
    double h1_error = 0.0;
    for (int cell = 0; cell < size; ++cell) {
        const int n_q_points = degree + 1;
        const int n_dofs = degree + 1;

        // evaluate the transform jacobian
        double jac_transform = dx / 2.0;

        // evaluate the shape functions and quadrature points
        shape_t shape_funcs[MAX_DEGREE][MAX_DEGREE];
        quadrature_t quadratures[MAX_DEGREE];

        for (int qp=0; qp<n_q_points; ++qp) {
            quadrature_t q = quadrature(qp, degree + 1);

            double u_val_num = 0.0;
            double u_grad_num = 0.0;
            for (int i=0; i<n_dofs; ++i) {
                shape_t shape = basis_eval(basis, q.point, i);
                shape.grad /= jac_transform;
                u_val_num += shape.value * u[cell * dofs_per_cell + i];
                u_grad_num += shape.grad * u[cell * dofs_per_cell + i];
            }

            double x0 = cell * dx;
            double x1 = x0 + dx;
            double x = (x0 + x1) * 0.5 + q.point * (x1 - x0) * 0.5;
            double dxq = q.weight * jac_transform;

            // exact solution
            #if defined(PROBLEM_POISSON)
                // poisson
                double freq = 6.28318530718;
                double u_val_exact = sin(freq * x);
                double u_grad_exact = freq * cos(freq * x);

            #elif defined(PROBLEM_CONVECT)
                // convection
                double u_val_exact = 2.0 * x * x / (x + 1.0) - x*x*x / 2.0;
                double u_grad_exact = x * (-3.0 * x*x*x - 6.0 * x*x + x + 8.0) / (2.0 * (x + 1.0) * (x + 1.0));
            #else
                #error "Define PROBLEM_POISSON or PROBLEM_CONVECT"
            #endif

            l2_error += (u_val_num - u_val_exact) * (u_val_num - u_val_exact) * dxq;
            h1_error += (u_grad_num - u_grad_exact) * (u_grad_num - u_grad_exact) * dxq;
        }
    }
    l2_error = sqrt(l2_error);
    h1_error = sqrt(h1_error);
    printf("- l2 error = %.5e\n", l2_error);
    printf("- h1 error = %.5e\n", h1_error);

    FILE* fe = fopen("plot/error.txt", "a");
    fprintf(fe, "%.5e, %.5e, %.5e\n", 1.0/dx, l2_error, h1_error);
    fclose(fe);


    printf("- Writing solution to files in plot/ dir\n");
    FILE* f = fopen("plot/data.txt", "w");
    FILE* fp = fopen("plot/edges.txt", "w");

    int n_p_per_elem = 10;
    if (size < 10) n_p_per_elem = 100 / size;
    for (int cell = 0; cell < size; ++cell) {
        double t[64];
        double y[64];
        for (int i=0; i<n_p_per_elem; ++i) {
            t[i] = -1.0 + 2.0 * ((double)i)/((double)(n_p_per_elem - 1));
            y[i] = 0.0;
            for (int j=0; j<dofs_per_cell; ++j) {
                y[i] += basis_eval(basis, t[i], j).value * u[cell * dofs_per_cell + j];
            }
        }
        double x0 = cell * dx;
        double x1 = x0 + dx;
        for (int i=0; i<n_p_per_elem; ++i)
            fprintf(f, "%lf, %lf\n", (x0 + x1) * 0.5 + t[i] * (x1 - x0) * 0.5, y[i]);
        
            fprintf(f, "%lf, nan\n", x1);
        
        fprintf(fp, "%lf %lf\n", x0, y[0]);
        fprintf(fp, "%lf %lf\n", x1, y[n_p_per_elem - 1]);
    }

    fclose(f);
    fclose(fp);

    // free linear system
    sparse_matrix_free(full_lhs);
    free(full_rhs);
    free(u);

    return 0;
}


// include the assembly functions

#if defined(PROBLEM_POISSON)
    #include "eqns/poisson.h"
#elif defined(PROBLEM_CONVECT)
    #include "eqns/convect.h"
#else
    #error "Define PROBLEM_POISSON or PROBLEM_CONVECT"
#endif

