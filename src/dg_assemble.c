#include "dg_assemble.h"
#include <stdio.h>
#include <stdlib.h>


#define MAX_DEGREE 12



void assemble_local_system(
    double** local_lhs,
    double* local_rhs,
    basis_t basis, // basis function
    double dx,  // cell width
    int* boundaries,  // boundary of (left face, right_face)
    int id,
    // assembly functions
    assembly_internal_lhs_t internal_lhs,
    assembly_internal_rhs_t internal_rhs,
    assembly_flux_t flux,
    assembly_boundary_lhs_t* boundary_lhs,
    assembly_boundary_rhs_t* boundary_rhs,
    void* user_data
) {
    int degree = basis.degree;

    if (degree > MAX_DEGREE) {
        fprintf(stderr, "error in assemble_local_system, degree (%d) > MAX_DEGREE (%d)", degree, MAX_DEGREE);
        exit(1);
    }

    const int n_q_points = degree + 1;
    const int n_dofs = degree + 1;

    // evaluate the transform jacobian
    double jac_transform = dx / 2.0;

    // evaluate the shape functions and quadrature points
    shape_t shape_funcs[MAX_DEGREE][MAX_DEGREE];
    quadrature_t quadratures[MAX_DEGREE];

    double x0 = ((double)id) * dx;
    double x1 = x0 + dx;

    for (int q=0; q<n_q_points; ++q) {
        quadratures[q] = quadrature(q, degree);
        for (int i=0; i<n_dofs; ++i) {
            shape_funcs[q][i] = basis_eval(basis, quadratures[q].point, i);
            shape_funcs[q][i].grad /= jac_transform;
        }
    }

    // compute the dxq integration weights
    double dxq[MAX_DEGREE];
    for (int q=0; q<n_q_points; ++q) {
        dxq[q] = quadratures[q].weight * jac_transform;
    }

    // initialize the local system
    for (int i = 0; i < n_dofs; ++i) {
        local_rhs[i] = 0.0;
        for (int j = 0; j < n_dofs; ++j) {
            local_lhs[0][i*n_dofs + j] = 0.0;
            if (local_lhs[1] != NULL) local_lhs[1][i*n_dofs + j] = 0.0;
            if (local_lhs[2] != NULL) local_lhs[2][i*n_dofs + j] = 0.0;
        }
    }

    // local system (self - self)
    for (int q = 0; q < n_q_points; ++q) {
        assembly_data_t data;
        data.q = quadratures[q].point;
        data.x = (x0 + x1) * 0.5 + data.q * (x1 - x0) * 0.5;

        for (int i = 0; i < n_dofs; ++i) {
            shape_t v = shape_funcs[q][i];
            for (int j = 0; j < n_dofs; ++j) {
                shape_t u = shape_funcs[q][j];

                local_lhs[0][i*n_dofs + j] += internal_lhs(u, v, dxq[q], data, user_data);
            }
            local_rhs[i] += internal_rhs(v, dxq[q], data, user_data);
        }
    }

    // face terms
    shape_t shape_0[2][MAX_DEGREE];
    shape_t shape_1[2][MAX_DEGREE];
    for (int i=0; i<n_dofs; ++i) {
        // left face
        shape_0[0][i] = basis_eval(basis, -1.0, i);
        shape_1[0][i] = basis_eval(basis, 1.0, i);

        shape_0[0][i].grad /= jac_transform;
        shape_1[0][i].grad /= jac_transform;

        // right face
        shape_0[1][i] = basis_eval(basis, 1.0, i);
        shape_1[1][i] = basis_eval(basis, -1.0, i);

        shape_0[1][i].grad /= jac_transform;
        shape_1[1][i].grad /= jac_transform;
    }


    double n[] = {-1.0, 1.0};
    double ds[] = {1.0, 1.0};
    assembly_data_t data[2];
    data[0].q = -1.0;
    data[0].x = x0;
    data[1].q = 1.0;
    data[1].x = x1;

    
    for (int f = 0; f < 2; ++f) {
        for (int i = 0; i < n_dofs; ++i) {
            shape_t v0 = shape_0[f][i];
            
            //shape_t v1 = shape_1[f][i];
            // always zero
            shape_t v1;
            v1.value = 0.0;
            v1.grad = 0.0;

            for (int j = 0; j < n_dofs; ++j) {
                shape_t u0, u1;

                if (boundaries[f] >= 0) {
                    // boundary face
                    u0 = shape_0[f][j];

                    local_lhs[0][i*n_dofs + j] += boundary_lhs[boundaries[f]](
                        u0, v0,
                        n[f], ds[f], dx, 
                        data[f],
                        user_data
                    );
                } else {
                    // internal face

                    // diagonal term
                    u0 = shape_0[f][j];
                    u1.value = 0.0;
                    u1.grad = 0.0;
                    local_lhs[0][i*n_dofs + j] += flux(
                        u0, v0, u1, v1,
                        n[f], ds[f], dx, 
                        data[f],
                        user_data
                    );

                    // off diagonal term
                    u1 = shape_1[f][j];
                    u0.value = 0.0;
                    u0.grad = 0.0;
                    local_lhs[f + 1][i*n_dofs + j] += flux(
                        u0, v0, u1, v1,
                        n[f], ds[f], dx, 
                        data[f],
                        user_data
                    );
                }
            }
            
            // rhs term
            if (boundaries[f] >= 0) {
                // this is a boundary face
                local_rhs[i] += boundary_rhs[boundaries[f]](
                    v0,
                    n[f], ds[f], dx, 
                    data[f],
                    user_data
                );
            } else {
                // this is an internal face
                // rhs term is null here
            }
        }
    }

    // done with assembly!
}
