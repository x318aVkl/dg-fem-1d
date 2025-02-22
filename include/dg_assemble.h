

#ifndef dg_assemble_h
#define dg_assemble_h

#include "dg_basis.h"


typedef struct assembly_data_t {
    double q;
    double x;
} assembly_data_t;

typedef double (*assembly_internal_lhs_t)(shape_t, shape_t, double, assembly_data_t, void* user_data);
typedef double (*assembly_internal_rhs_t)(shape_t, double, assembly_data_t, void* user_data);
typedef double (*assembly_flux_t)(shape_t, shape_t, shape_t, shape_t, double, double, double, assembly_data_t, void* user_data);
typedef double (*assembly_boundary_lhs_t)(shape_t, shape_t, double, double, double, assembly_data_t, void* user_data);
typedef double (*assembly_boundary_rhs_t)(shape_t, double, double, double, assembly_data_t, void* user_data);



// assembles the linear system for one cell
// local system must be of size lhs(degree + 1, degree + 1), rhs(degree+1)
// local matrices must point to three matrices
// local_matrices[0] is diagonal, local_matrices[1] is left cell, local_matrices[2] is right cell
void assemble_local_system(
    double** local_lhs,
    double* local_rhs,
    basis_t basis,  // basis function
    double dx,  // cell width
    int* boundaries,  // boundary of (left face, right_face)
    int id, // cell id
    // assembly functions
    assembly_internal_lhs_t internal_lhs,
    assembly_internal_rhs_t internal_rhs,
    assembly_flux_t flux,
    assembly_boundary_lhs_t* boundary_lhs,
    assembly_boundary_rhs_t* boundary_rhs,
    void* user_data
);



#endif
