// The file unimodular_matrices.h specifies functions to generate
// random unimodular matrices.

#ifndef __unimodular_matrices_h__
#define __unimodular_matrices_h__

void random_row_and_column ( int dim, int *row, int *col );
/*
 * DESCRIPTION :
 *   Returns a pair of two different indices (row, col),
 *   in the range 0..dim-1. */

int nonzero_integer ( int size, bool poscff=false );
/*
 * DESCRIPTION :
 *   Returns a number in the range -size .. + size,
 *   if poscff then the returned number is positive. */

void unimodular_multiplication
 ( int dim, int row, int col, int nbr, int **mat );
/*
 * DESCRIPTION :
 *   Updates the matrix in mat with the product of the unimodular
 *   matrix of dimension dim defined by the number nbr at (row, col). */

int maximum ( int dim, int **mat );
/*
 * DESCRIPTION :
 *   Returns the largest element in the matrix mat of dimension dim. */

void make_unimodular_matrix
 ( int dim, int size, bool poscff, int itr, int **mat, bool verbose );
/*
 * DESCRIPTION :
 *   Generates a unimodular matrix of dimension dim,
 *   with number of the given size, all positive if poscff,
 *   using itr iterations of multiplications with a unimodular matrix.
 *   On input is memory allocated in the matrix mat.
 *   The result is returned in mat. */

void write_exponent_matrix ( int dim, int **mat );
/*
 * DESCRIPTION :
 *   Writes the matrix of dimension dim defined by mat. */

void read_exponent_matrix ( int dim, int **mat );
/*
 * DESCRIPTION :
 *   Prompts the user for a matrix of dimension dim,
 *   which is returned in the matrix mat. */

void lower_triangular_unit ( int dim, int **mat );
/*
 * DESCRIPTION :
 *   Fills up the matrix of dim rows and dim colums
 *   with zeros above the diagonal and ones on and below the diagonal. */

void two_variable_monomials ( int dim, int **mat );
/*
 * DESCRIPTION :
 *   The matrix on return corresponds to a decoupled system
 *   where each monomial has exactly two participating variables,
 *   except the first one in case the dimension is odd. */

void make_monomial_system
 ( int dim, int size, int posvals, int nbritr,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   int **rowsA, int vrblvl );
/*
 * DESCRIPTION :
 *   Generates a monomial system with a unimodular exponent matrix.
 *
 * ON ENTRY :
 *   dim       the number of monoimials;
 *   size      size of the random exponents;
 *   posvals   if positive, then all exponents are positive;
 *   nbritr    number of unimodular multiplications,
 *             if zero, then the user is prompted for a matrix;
 *   nvr       space for dim many integers;
 *   idx       space for dim many pointers to integers;
 *   exp       space for dim many pointers to integers;
 *   nbrfac    space for dim many integers;
 *   expfac    space for dim many pointers to integers;
 *   rowsA     space for dim many pointers to integers;
 *   vrblv     is the verbose level, 0 if silent.
 *
 * ON RETURN :
 *   nvr       nvr[i] is the number of variables in the i-th monomial;
 *   idx       idx[i] is an array of size nvr[i],
 *             with the indices of the variables in the i-th monomial;
 *   exp       exp[i] is an array of size nvr[i],
 *             with the exponents of the variables in the i-th monomial;
 *   nbrfac    nbrfac[i] has the number of exponents that are
 *             larger than one in the i-th monomial;
 *   expfac    expfac[i] are the exponents in the i-th polynomial
 *             that are larger than one, minus one in the factor,
 *             if exp[i][k] > 1, then expfac[i][k] = exp[i][k] - 1;
 *   rowsA     contains the dim rows of the exponent matrix. */

void prompt_dimensions
 ( int *dim, int *deg, int *size, int *posvals, int *vrblvl, int *nbritr,
   int *nbsteps );
/*
 * DESCRIPTION :
 *   Prompts for the dimensions.
 *
 * ON RETURN :
 *   dim       the dimension is the number of monomials
 *             and the maximum number of variables in each monomial;
 *   deg       degree of the series;
 *   size      size of the numbers;
 *   posvals   positive exponents (1 if yes);
 *   vrblvl    verbose level (0 if silent);
 *   nbritr    number of unimodular multiplications in the making of
 *             the exponent matrix, otherwise, for nonpositive values:
 *             if -3, then cyclic n-roots,
 *             if -2, then decoupled two variable monomials,
 *             if -1, then lower triangular matrix of ones,
 *             if 0, then user input;
 *   nbsteps   the number of Newton steps. */

void copy_integer_matrix ( int dim, int **matfrom, int **matinto );
/*
 * DESCRIPTION :
 *   Copies the matrix matfrom into the matrix matinto.
 *
 * ON ENTRY :
 *   dim       the dimension of the matrices matfrom and matinto;
 *   matfrom   originating matrix;
 *   matinto   space allocated for a matrix of dimensiom dim.
 *
 * ON RETURN :
 *   matinto   the copy of the matrix matfrom. */

void matrix_matrix_multiply ( int dim, int **A, int **B, int **C );
/*
 * DESCRIPTION :
 *   Multiplies A with B and stores the result in C.
 *   All matrices are square and of dimension dim. */

int posgcd ( int a, int b, int *ka, int *lb, int vrblvl );
/*
 * DESCRIPTION :
 *   Returns the greatest common divisor d of a and b
 *   and the cofactors, i.e.: d = ka*a + lb*b,
 *   for positive inputs a and b.
 *
 * ON ENTRY :
 *   a,b       two positive integer numbers;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   ka        cofactor to multiply a with,
 *   lb        cofactor to multiply b with. */

int gcd ( int a, int b, int *ka, int *lb, int vrblvl );
/*
 * DESCRIPTION :
 *   Returns the greatest common divisor d of a and b
 *   and the cofactors, i.e.: d = ka*a + lb*b.
 *
 * ON ENTRY :
 *   a,b       two integer numbers;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   ka        cofactor to multiply a with,
 *   lb        cofactor to multiply b with. */

int lower_triangulate ( int dim, int **mat, int **uni, int vrblvl );
/*
 * DESCRIPTION :
 *   Transforms the matrix into a lower triangular matrix,
 *   via unimodular coordinate transformations.
 *   Returns zero if the matrix is regular,
 *   -1 if a zero column was computed and the matrix is singular.
 *
 * ON ENTRY :
 *   dim       the dimension of the matrices, number of rows and columns;
 *   mat       input matrix of dimension dim;
 *   uni       space allocated for a matrix of dimension dim;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   mat       is an upper triangular matrix;
 *   uni       the unimodular coordinate transformation matrix:
 *             uni multiplied with the original matrix equals mat. */

int exponent_forward_substitution ( int dim, int **mat, int *expsol );
/*
 * DESCRIPTION :
 *   If all right hand side elements of a binomial system equal the same c,
 *   then the components of the solution are powers of that constant c.
 *   Those powers are computed for an lower triangular exponent matrix.
 *   Returns -1 if the diagonal of the matrix has a zero, else returns 0.
 *
 * ON ENTRY :
 *   dim       dimension of the exponent matrix;
 *   mat       upper triangular matrix for forward substitution;
 *   expsol    space for dim integers.
 *
 * ON RETURN :
 *   expsol    exponents of the constant in the solution vector. */

int exponent_unimodular_transformation ( int dim, int **uni, int *expsol );
/*
 * DESCRIPTION :
 *   Following the exponent forward substitution, the exponents are
 *   subjected to the unimodular coordinate transformation used to
 *   bring the exponent matrix into lower triangular form.
 *
 * ON ENTRY :
 *   dim       dimension of the matrix and the exponent vector;
 *   uni       unimodular coordinate transformation;
 *   expsol    are the current values of the exponents of the solution.
 *
 * ON RETURN :
 *   expsol    exponents of the constant c in the right hand side
 *             of a binomial system, to define the solution,
 *             obtained via multiplying uni by expsol. */

int exponents_check
 ( int dim, int **rowsA, int *expsol, int vrblvl );
/*
 * DESCRIPTION :
 *   Writes the exponents of the solution, returned in expsol,
 *   for the matrix A stored in rowsA of dimension dim.
 *   The vector in expsol has dimension dim.
 *   The return value is zero if the matrix in rowsA is nonsingular,
 *   otherwise the return value is -1.
 *   If vrblvl > 1, then the matrix is written as well. */

void row_sums ( int dim, int **rowsA, int *sums );
/*
 * DESCRIPTION :
 *   Computes the sums of the rows in the matrix.
 *
 * ON ENTRY :
 *   dim      number of rows in the matrix and length of sums;
 *   rowsA    matrix of dimension dim stored row wise;
 *   sums     space for dim integers.
 *
 * ON RETURN :
 *   sums     sums[k] has the sum of the k-th row of rowsA. */

#endif
