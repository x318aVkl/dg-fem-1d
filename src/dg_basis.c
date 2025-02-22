

#include "dg_basis.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>



// shape function defined on the element
// -1 ----- t ------ 1
// t position, i id, d degree + 1

// legendre polynomials
double legendre_shape_value(double t, int i) {
    #include "basis/legendre_shape.h"
    fprintf(stderr, "shape_function(): invalid id %d\n", i);
    exit(1);
    return 0;
}

double legendre_shape_grad(double t, int i) {
    #include "basis/legendre_grad.h"
    fprintf(stderr, "shape_grad(): invalid id %d\n", i);
    exit(1);
    return 0;
}



// lagrange polynomials
double lagrange_shape_value(double t, int i, int d) {
    #include "basis/lagrange_shape.h"
    fprintf(stderr, "shape_function(): invalid id %d or degree %d\n", i, d);
    exit(1);
    return 0;
}

double lagrange_shape_grad(double t, int i, int d) {
    #include "basis/lagrange_grad.h"
    fprintf(stderr, "shape_grad(): invalid id %d or degree %d", i, d);
    exit(1);
    return 0;
}



double quadrature_point(int q, int d) {
    #include "basis/quadrature_point.h"
    printf("ERROR in quadrature_point(): invalid id %d or degree %d\n", q, d);
    exit(1);
    return 0;
}

double quadrature_weight(int q, int d) {
    #include "basis/quadrature_weight.h"
    printf("ERROR in quadrature_point(): invalid id %d or degree %d\n", q, d);
    exit(1);
    return 0;
}


quadrature_t quadrature(
    int q_point,
    int degree
) {
    quadrature_t q;

    q.point = quadrature_point(q_point, degree);
    q.weight = quadrature_weight(q_point, degree);

    return q;
}




shape_t legendre_shape_function(
    double t,
    int dof,
    int _degree
) {
    shape_t s;

    s.value = legendre_shape_value(t, dof);
    s.grad = legendre_shape_grad(t, dof);

    return s;
}

shape_t lagrange_shape_function(
    double t,
    int dof,
    int degree
) {
    shape_t s;

    s.value = lagrange_shape_value(t, dof, degree);
    s.grad = lagrange_shape_grad(t, dof, degree);

    return s;
}




basis_t legendre_basis(int degree) {
    basis_t basis;

    basis.degree = degree;
    basis.shape = legendre_shape_function;

    return basis;
}


basis_t lagrange_basis(int degree) {
    basis_t basis;

    basis.degree = degree;
    basis.shape = lagrange_shape_function;

    return basis;
}




shape_t basis_eval(basis_t basis, double t, int i) {
    return basis.shape(t, i, basis.degree);
}


