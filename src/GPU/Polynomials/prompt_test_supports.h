/* The file prompt_test_supports.h specifies functions to generate
 * the supports of test polynomials. */

void prompt_testpoly_dimensions ( int *seed, int *dim, int *nva, int *nbr );
/*
 * DESCRIPTION :
 *   Prompts for the dimensions of the test polynomial.
 *
 * ON RETURN :
 *   seed      seed used for the random number generator;
 *   dim       the dimension is the total number of variables;
 *   nva       number of variables per monomial, 0 for a random polynomial;
 *   nbr       number of terms in the polynomial. */

void make_test_supports ( int dim, int nva, int nbr, int *nvr, int **idx );
/*
 * DESCRIPTION :
 *   Generates the supports of a test polynomial.
 *
 * ON ENTRY :
 *   dim       total number of variables;
 *   nva       number of variables per monomial;
 *   nbr       number of terms in the polynomial;
 *   nvr       space for nbr integers;
 *   idx       space for nbr pointers to integers.
 *
 * ON RETURN :
 *   nvr       number of variables in each monomial;
 *   idx       indices of the variables in each monomial. */
