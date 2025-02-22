
#ifndef spmat_h
#define spmat_h



// a simple row major sparse matrix type
// with a fixed number of entries per rows
typedef struct {
    // actual values
    int* columns;
    double* values;

    // size info
    int n_rows;
    int n_per_row;
} sparse_matrix_t;



// create a sparse matrix
sparse_matrix_t sparse_matrix_create(
    int nrows,
    int n_per_row
);


// clean up memory
void sparse_matrix_free(
    sparse_matrix_t spmat
);


// index sparse matrix
double* sparse_matrix_index(
    sparse_matrix_t spmat,
    int i,
    int j
);


// solve linear problem using gauss siedel
// with ssor, omega parameter, if unsure set to 1
// returns the number of iterations
// if solving a non diagonaly dominant problem, use an omega value smaller than 1, ex. 0.5
int gauss_siedel(
    sparse_matrix_t lhs,
    const double* rhs,
    double* x,
    double omega,
    double tol,
    int verbose
);





#endif

