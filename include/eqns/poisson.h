

// assembly functions for the poisson equation with zero dirichlet bcs


// assembly functions
double interior_lhs(
    shape_t u,
    shape_t v,
    double dx,
    assembly_data_t data,
    void* user_data
) {
    return v.grad * u.grad * dx;
}


double interior_rhs(
    shape_t v,
    double dx,
    assembly_data_t data,
    void* user_data
) {
    // manufactured solution
    // u = sin(2 * pi * x)
    // u = sin(f * pi)
    // u'' = - f*f * sin(f * x)
    double f = 6.28318530718;
    return f*f * sin(f * data.x) * v.value * dx;
}



//  |   other cell(1)     <-n-|   this cell(0)  |
double numerical_flux(
    shape_t u0,
    shape_t v0,
    shape_t u1,
    shape_t v1,
    double n,
    double ds,
    double h, 
    assembly_data_t data,
    void* user_data
) {
    int degree = *((int*)user_data);
    double delta = -1.0;
    double order = ((double)degree);
    double alpha_stab = order * (order + 1) / h;
    if (degree == 0) {
        alpha_stab = 1.0 / h;
    }

    return (
          (v1.value - v0.value) * 0.5 * (u1.grad*n + u0.grad*n)
        - delta * (v1.grad*n + v0.grad*n) * 0.5 * (u1.value - u0.value)
        + alpha_stab * (v1.value - v0.value) * (u1.value - u0.value)
    ) * ds;
}




double boundary_term_lhs(
    shape_t u,
    shape_t v,
    double n,
    double ds,
    double h, 
    assembly_data_t data,
    void* user_data
) {
    int degree = *((int*)user_data);
    double delta = -1.0;
    double order = ((double)degree);
    double alpha_stab = order * (order + 1) / h;
    if (degree == 0) {
        alpha_stab = 1.0 / h;
    }

    return (
        ( - v.value ) * 0.5 * ( 2.0 * u.grad*n )
        - delta * ( v.grad*n ) * 0.5 * ( - u.value )
        + alpha_stab * ( - v.value) * ( - u.value)
    ) * ds;
}


double boundary_term_rhs(
    shape_t v,
    double n,
    double ds,
    double h, 
    assembly_data_t data,
    void* user_data
) {
    int degree = *((int*)user_data);

    double dirichlet = 0.0;

    double delta = -1.0;
    double order = ((double)degree);
    double alpha_stab = order * (order + 1) / h;
    if (degree == 0) {
        alpha_stab = 1.0 / h;
    }
    return - (
        - delta * ( v.grad*n ) * 0.5 * ( dirichlet )
        + alpha_stab * ( - v.value) * ( dirichlet )
    ) * 1.0 * ds;
}




