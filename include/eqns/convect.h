

// assembly functions for the convection/diffusion equation with zero dirichlet bcs


// assembly functions
double interior_lhs(
    shape_t u,
    shape_t v,
    double dx, 
    assembly_data_t data,
    void* user_data
) {
    return - v.grad * u.value * dx;
}


double interior_rhs(
    shape_t v,
    double dx, 
    assembly_data_t data,
    void* user_data
) {
    double x = data.x;
    double f = x * (-3.0 * x*x*x - 6.0 * x*x + x + 8.0) / (2.0 * (x + 1.0) * (x + 1.0));
    return f * v.value * dx;
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
    // simple upwind flux

    if (n > 0.0) {
        return v0.value * u0.value * n * ds;
    } else {
        return v0.value * u1.value * n * ds;
    }
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
    if (n > 0.0) {
        return v.value * u.value * n * ds;
    } else {
        return 0.0; //- v.value * u.value * n * ds;
    }
}


double boundary_term_rhs(
    shape_t v,
    double n,
    double ds,
    double h, 
    assembly_data_t data,
    void* user_data
) {
    if (n > 0.0) {
        // zero flux
        return v.value * 0.0 * n * ds;
    } else {
        return v.value * 0.0 * n * ds;
    }
}




