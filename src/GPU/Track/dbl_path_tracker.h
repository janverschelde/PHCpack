// The file dbl_path_tracker.h specifies a path tracker
// on series in double precision on real numbers.

#ifndef __dbl_path_tracker_h__
#define __dbl_path_tracker_h__

int dbl_run_newton
 ( int szt, int nbt, int dim, int deg, int nbrcol, int nbsteps,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac,
   double **mbrhs, double ***cff, double **acc,
   double **input_h, double **input_d, double ***output_h, double ***output_d,
   double **funval_h, double **funval_d,
   double ***jacval_h, double ***jacval_d, double **rhs_h, double **rhs_d,
   double **urhs_h, double **urhs_d, double **sol_h, double **sol_d,
   double **Q_h, double **Q_d, double **R_h, double **R_d,
   double *workvec, double **resvec, double *resmax,
   int vrblvl, int mode );
/*
 * DESCRIPTION :
 *   Runs Newton's method on power series on real data in double precision.
 *
 * REQUIRED : szt*nbt = dim for GPU computing.
 *
 * ON ENTRY :
 *   szt       size of each tile and block;
 *   nbt       number of tiles and number of blocks;
 *   dim       number of monomials;
 *   deg       degree of the power series;
 *   nbrcol    is the number of columns, if 1, then the system is monomial,
 *             otherwise nbrcol columns are given on input;
 *   nbsteps   maximum number of Newton steps;
 *   nvr       nvr[i][j] is the number of variables in the j-th monomial
 *             of the i-th column;
 *   idx       idx[i][j] are the indices of the variables in monomial j
 *             of the i-th column;
 *   exp       exp[i] are the exponents of the variables in monomial i;
 *   nbrfac    nbrfac[i] are the number of exponents > 1 in monomial i;
 *   expfac    expfac[i] are the exponents in the i-th polynomial
 *             that are larger than one, minus one in the factor,
 *             if exp[i][k] > 1, then expfac[i][k] = exp[i][k] - 1;
 *   mbrhs     right hand side vector of series;
 *   cff       coefficients of the monomials in each column;
 *   acc       space to accumulate dim+1 power series of degree deg;
 *   input_h   coefficients of the power series of degree deg,
 *             for dim variables, computed on host;
 *   input_d   space for power series computed on device;
 *   output_h  space for the evaluated and differentiated monomials,
 *             computed on the host;
 *   output_d  space for the evaluated and differentiated monomials,
 *             computed on the device;
 *   funval_h  space for the evaluated power series computed by host;
 *   funval_d  space for the evaluated power series computed by device;
 *   jacval_h  space for deg+1 matrices of dimension dim on host;
 *   jacval_d  space for deg+1 matrices of dimension dim on device;
 *   rhs_h     space for deg+1 vectors of dimension dim on host;
 *   rhs_d     space for deg+1 vectors of dimension dim on device;
 *   urhs_h    space for updated right hand side vectors computed by host;
 *   urhs_d    space for updated right hand side vectors computed by device; 
 *   sol_h     space for deg+1 vectors of dimension dim;
 *   sol_d     space for deg+1 vectors of dimension dim;
 *   Q_h       space allocated for the Q computed by the host;
 *   Q_d       space allocated for the Q computed by the device;
 *   R_h       space allocated for the R computed by the host;
 *   R_d       space allocated for the R computed by the device;
 *   wrkvec    work space allocated for a vector of dimension dim;
 *   resvec    space for deg+1 vectors of dimension dim;
 *   vrblvl    is the verbose level;
 *   mode      execution mode, 0 (GPU only), 1 (CPU only) or 2 (GPU+CPU).
 *
 * ON RETURN :
 *   input_h   power series computed on host (depending on mode);
 *   input_d   power series computed on device (depending on mode);
 *   output_h  evaluated power series computed on host (depending on mode);
 *   output_d  evaluated power series computed on device (depending on mode);
 *   funval_h  collects the output[i][dim], the evaluated series on host;
 *   funval_d  collects the output[i][dim], the evaluated series on device;
 *   jacval_h  a series with matrices as coefficients, computed by host,
 *             the leading coefficient is the Jacobian matrix;
 *   jacval_d  a series with matrices as coefficients, computed by device,
 *             the leading coefficient is the Jacobian matrix;
 *   rhs_h     the linearized right hand side are the function values
 *             subtracted by 1 and added by t, computed by host;
 *   rhs_d     the linearized right hand side are the function values
 *             subtracted by 1 and added by t, computed by device;
 *   urhs_h    right hand side vector updated by the host;
 *   urhs_d    right hand side vector updated by the device;
 *   sol_h     solution computed by the host;
 *   sol_d     solution computed by the device;
 *   Q_h       Q of the QR factorization computed by the host;
 *   Q_d       Q of the QR factorization computed by the device;
 *   R_h       R of the QR factorization computed by the host;
 *   R_d       R of the QR factorization computed by the device;
 *   wrkmat    has a copy of the Jacobian matrix;
 *   resvec    residual vectors;
 *   resmax    the maximum element of the residual vectors. */

int test_dbl_real_track
 ( int szt, int nbt, int dim, int deg, int nbrcol,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac, int **rowsA,
   int nbsteps, int mode, int vrblvl );
/*
 * DESCRIPTION :
 *   Tracks a path on a system with real double arithmetic.
 *
 * ON ENTRY :
 *   szt       size of each tile and block;
 *   nbt       number of tiles and number of blocks;
 *   dim       number of monomials;
 *   deg       degree of the power series;
 *   nbrcol    number of columns, if 1, then the system is monomial,
 *             otherwise nbrcol columns are expected;
 *   nvr       nvr[i][j] is the number of variables in the j-th monomial
 *             of the i-th column;
 *   idx       idx[i][j] are the indices of the variables in monomial j
 *             of the i-th column;
 *   exp       exp[i] are the exponents of the variables in monomial i;
 *   nbrfac    nbrfac[i] are the number of exponents > 1 in monomial i;
 *   expfac    expfac[i] are the exponents in the i-th polynomial
 *             that are larger than one, minus one in the factor,
 *             if exp[i][k] > 1, then expfac[i][k] = exp[i][k] - 1;
 *   rowsA     rows of the exponents of the dim monomials;
 *   nbsteps   the number of Newton steps;
 *   mode      the mode of execution, 0 for GPU only, 1 for CPU only,
 *             2 for CPU and GPU;
 *   vrblvl    is the verbose level. */

#endif
