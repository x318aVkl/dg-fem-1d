
#ifndef dg_basis_h
#define dg_basis_h



// shape function evaluation
typedef struct {
    double value;
    double grad;
} shape_t;


typedef struct {
    double point;
    double weight;
} quadrature_t;



typedef shape_t (*shape_function)(double, int, int);


typedef struct {
    shape_function shape;
    int degree;
} basis_t;



//shape_t shape_function(
//    double t,
//    int dof
//);


basis_t legendre_basis(int degree);
basis_t lagrange_basis(int degree);

shape_t basis_eval(basis_t basis, double t, int i);


quadrature_t quadrature(
    int q_point,
    int degree
);





#endif

