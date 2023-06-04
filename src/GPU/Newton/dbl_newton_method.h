// The file dbl_newton_method.h specifies Newton's method
// on series in double precision on real numbers.

#ifndef __dbl_newton_method_h__
#define __dbl_newton_method_h__

int dbl_errors_funjacrhs
 ( int dim, int deg, double **funval_h, double **funval_d,
   double ***jacval_h, double ***jacval_d,
   double **rhs_h, double **rhs_d, int vrblvl );
/*
 * DESCRIPTION :
 *   Compares results on host with those on the device
 *   and returns the number of errors larger than 1.0e-8.
 *
 * ON ENTRY :
 *   dim       number of equations and variables;
 *   deg       degree at which the series are truncated;
 *   funval_h  function values on the host;
 *   funval_d  function values on the device;
 *   jacval_h  values of the Jacobian matrix on the host;
 *   jacval_d  values of the Jacobian matrix on the device;
 *   rhs_h     right hand side vector on the host;
 *   rhs_d     right hand side vector on the device;
 *   vrblvl    verbose level. */

int dbl_errors_inurhsQRsol
 ( int dim, int deg, double **input_h, double **input_d,
   double **Q_h, double **Q_d, double **R_h, double **R_d,
   double **urhs_h, double **urhs_d, double **sol_h, double **sol_d,
   int vrblvl );
/*
 * DESCRIPTION :
 *   Compares results on host with those on the device
 *   and returns the number of errors larger than 1.0e-8.
 *
 * ON ENTRY :
 *   dim       number of equations and variables;
 *   deg       degree at which the series are truncated;
 *   input_h   updated series on the host;
 *   input_d   updated series on the device;
 *   Q_h       Q on the host;
 *   Q_d       Q on the device;
 *   R_h       R on the host;
 *   R_d       R on the device;
 *   urhs_h    updated right hand side on the host;
 *   urhs_d    updated right hand side on the device;
 *   sol_h     update to the solution on the host,
 *   sol_d     update to the solution on the device;
 *   vrblvl    is the verbose level. */

int dbl_update_newton_qrstep
 ( int szt, int nbt, int dim, int deg,
   int *tailidx_h, int *tailidx_d, double **input_h, double **input_d,
   double **funval_h, double **funval_d,
   double ***jacval_h, double ***jacval_d, double **rhs_h, double **rhs_d,
   double **urhs_h, double **urhs_d, double **sol_h, double **sol_d,
   double **Q_h, double **Q_d, double **R_h, double **R_d,
   double *workvec, double **resvec, double *resmax,
   bool *zeroQ_h, bool *noqr_h, bool *zeroQ_d, bool *noqr_d,
   int *upidx_h, int *bsidx_h, int *upidx_d, int *bsidx_d,
   double *totqrlapsedms, double *totqtblapsedms,
   double *totbslapsedms, double *totupdlapsedms, double *totreslapsedms,
   int vrblvl, int mode );
/*
 * DESCRIPTION :
 *   Given the function values and the matrix series of the Jacobian,
 *   computes the update to the solution series.
 *
 * REQUIRED : szt*nbt = dim for GPU acceleration, when mode is 2.
 *
 * ON ENTRY :
 *   szt       size of each tile and block;
 *   nbt       number of tiles and number of blocks;
 *   dim       number of monomials;
 *   deg       degree of the power series;
 *   tailidx_h is the start index of the update of the tail on the host;
 *   tailidx_d is the start index of the update of the tail on the device;
 *   input_h   coefficients of the power series of degree deg,
 *             for dim variables, computed on host;
 *   input_d   series computed on device;
 *   funval_h  evaluated power series computed by host;
 *   funval_d  evaluated power series computed by device;
 *   jacval_h  deg+1 matrices of dimension dim on host;
 *   jacval_d  deg+1 matrices of dimension dim on device;
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
 *   zeroQ_h   if true, then Q is zero and Q must be computed on host;
 *   noqr_h    flag if true, then no qr on host;
 *   zeroQ_d   if true, then Q is zero and Q must be computed on device;
 *   noqr_d    flag if true, then no qr on device;
 *   totqrlapsedms accumulates the milliseconds spent on the Householder QR;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs;
 *   totbslapsedms accumulates the milliseconds spent on back substitutions;
 *   totupdlapsedms accumulates the milliseconds spent on updates;
 *   totreslapsedms accumulates the milliseconds spent on residuals;
 *   vrblvl    is the verbose level;
 *   mode      execution mode, 0 (GPU only), 1 (CPU only) or 2 (GPU+CPU).
 *
 * ON RETURN :
 *   tailidx_h is the updated value for tailidx_h;
 *   tailidx_d is the updated value for tailidx_d;
 *   input_h   power series computed on host (depending on mode);
 *   input_d   power series computed on device (depending on mode);
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
 *   resmax    the maximum element of the residual vectors;
 *   zeroQ_h   false if Q was computed on host;
 *   noqr_h    updated flag if ||dx_0|| is zero for the first time on host;
 *   zeroQ_d   false if Q was computed on device;
 *   noqr_d    updated flag if ||dx_0|| is zero for the first time on device;
 *   upidx_h   counts the number of updates skipped by host;
 *   bsidx_h   counts the number of backsubstitutions skipped by host;
 *   upidx_d   counts the number of updates skipped by device;
 *   bsidx_d   counts the number of backsubstitutions skipped by device;
 *   totqrlapsedms accumulates the milliseconds spent on the Householder QR;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs;
 *   totbslapsedms accumulates the milliseconds spent on back substitutions;
 *   totupdlapsedms accumulates the milliseconds spent on updates;
 *   totreslapsedms accumulates the milliseconds spent on residuals. */

int dbl_column_newton_qrstep
 ( int szt, int nbt, int dim, int deg, int nbrcol,
   int *tailidx_h, int *tailidx_d,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac,
   double **mb, double dpr, double ***cff, double **acc,
   double **input_h, double **input_d, double ***output_h, double ***output_d,
   double **funval_h, double **funval_d,
   double ***jacval_h, double ***jacval_d, double **rhs_h, double **rhs_d,
   double **urhs_h, double **urhs_d, double **sol_h, double **sol_d,
   double **Q_h, double **Q_d, double **R_h, double **R_d,
   double *workvec, double **resvec, double *resmax,
   bool *zeroQ_h, bool *noqr_h, bool *zeroQ_d, bool *noqr_d, 
   int *upidx_h, int *bsidx_h, int *upidx_d, int *bsidx_d,
   double *totcnvlapsedms, double *totqrlapsedms, double *totqtblapsedms,
   double *totbslapsedms, double *totupdlapsedms, double *totreslapsedms,
   int vrblvl, int mode );
/*
 * DESCRIPTION :
 *   Does one step with Newton's method to update a power series,
 *   using QR factorization to solve linear systems, on real data,
 *   on one or more columns of monomials.
 *
 * REQUIRED : szt*nbt = dim for GPU acceleration, when mode is 2.
 *
 * ON ENTRY :
 *   szt       size of each tile and block;
 *   nbt       number of tiles and number of blocks;
 *   dim       number of monomials;
 *   deg       degree of the power series;
 *   nbrcol    is the number of columns, if 1, then the system is monomial,
 *             otherwise nbrcol columns are given on input;
 *   tailidx_h is the start index of the update of the tail on the host;
 *   tailidx_d is the start index of the update of the tail on the device;
 *   nvr       nvr[i][j] is the number of variables in the j-th monomial
 *             of the i-th column;
 *   idx       idx[i][j] are the indices of the variables in monomial j
 *             of the i-th column;
 *   exp       exp[i] are the exponents of the variables in monomial i;
 *   nbrfac    nbrfac[i] are the number of exponents > 1 in monomial i;
 *   expfac    expfac[i] are the exponents in the i-th polynomial
 *             that are larger than one, minus one in the factor,
 *             if exp[i][k] > 1, then expfac[i][k] = exp[i][k] - 1;
 *   mb        right hand side vector of series;
 *   dpr       damper multiplier for t, should be in (0.0, 1.0];
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
 *   zeroQ_h   if true, then Q is zero and Q must be computed on host;
 *   noqr_h    flag if true, then no qr on host;
 *   zeroQ_d   if true, then Q is zero and Q must be computed on device;
 *   noqr_d    flag if true, then no qr on device;
 *   totcnvlapsedms accumulates the milliseconds spent on the convolutions;
 *   totqrlapsedms accumulates the milliseconds spent on the Householder QR;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs;
 *   totbslapsedms accumulates the milliseconds spent on back substitutions;
 *   totupdlapsedms accumulates the milliseconds spent on updates;
 *   totreslapsedms accumulates the milliseconds spent on residuals;
 *   vrblvl    is the verbose level;
 *   mode      execution mode, 0 (GPU only), 1 (CPU only) or 2 (GPU+CPU).
 *
 * ON RETURN :
 *   tailidx_h is the updated value for tailidx_h;
 *   tailidx_d is the updated value for tailidx_d;
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
 *   resmax    the maximum element of the residual vectors;
 *   zeroQ_h   false if Q was computed on host;
 *   noqr_h    updated flag if ||dx_0|| is zero for the first time on host;
 *   zeroQ_d   false if Q was computed on device;
 *   noqr_d    updated flag if ||dx_0|| is zero for the first time on device;
 *   upidx_h   counts the number of updates skipped by host;
 *   bsidx_h   counts the number of backsubstitutions skipped by host;
 *   upidx_d   counts the number of updates skipped by device;
 *   bsidx_d   counts the number of backsubstitutions skipped by device;
 *   totcnvlapsedms accumulates the milliseconds spent on the convolutions;
 *   totqrlapsedms accumulates the milliseconds spent on the Householder QR;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs;
 *   totbslapsedms accumulates the milliseconds spent on back substitutions;
 *   totupdlapsedms accumulates the milliseconds spent on updates;
 *   totreslapsedms accumulates the milliseconds spent on residuals. */

int dbl_row_newton_qrstep
 ( int szt, int nbt, int dim, int deg, int *tailidx_h, int *tailidx_d,
   int *nbr, int **nvr, int ***idx,
   double **cst, double ***cff, double dpr,
   double **input_h, double **input_d, double ***output_h, double ***output_d,
   double **funval_h, double **funval_d,
   double ***jacval_h, double ***jacval_d, double **rhs_h, double **rhs_d,
   double **urhs_h, double **urhs_d, double **sol_h, double **sol_d,
   double **Q_h, double **Q_d, double **R_h, double **R_d,
   double *workvec, double **resvec, double *resmax,
   bool *zeroQ_h, bool *noqr_h, bool *zeroQ_d, bool *noqr_d,
   int *upidx_h, int *bsidx_h, int *upidx_d, int *bsidx_d,
   double *totcnvlapsedms, double *totqrlapsedms, double *totqtblapsedms,
   double *totbslapsedms, double *totupdlapsedms, double *totreslapsedms,
   int vrblvl, int mode );
/*
 * DESCRIPTION :
 *   Does one step with Newton's method to update a power series,
 *   using QR factorization to solve linear systems, on real data,
 *   on an indexed polynomial system.
 *
 * REQUIRED : szt*nbt = dim for GPU acceleration, when mode is 2.
 *
 * ON ENTRY :
 *   szt       size of each tile and block;
 *   nbt       number of tiles and number of blocks;
 *   dim       number of monomials;
 *   deg       degree of the power series;
 *   tailidx_h is the start index of the update of the tail on the host;
 *   tailidx_d is the start index of the update of the tail on the device;
 *   nbr       nbr[i] is the number of monomials in the i-th polynomial;
 *   nvr       nvr[i][j] is the number of variables in the j-th monomial
 *             of the i-th polynomial;
 *   idx       idx[i][j] are the indices of the variables in monomial j
 *             of the i-th polynomial;
 *   cst       vector of the constant coefficients of the system;
 *   cff       coefficients of the monomials in the system;
 *   dpr       damper multiplier for t, should be in (0.0, 1.0];
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
 *   zeroQ_h   if true, then Q is zero and Q must be computed on host;
 *   noqr_h    flag if true, then no qr on host;
 *   zeroQ_d   if true, then Q is zero and Q must be computed on device;
 *   noqr_d    flag if true, then no qr on device;
 *   totcnvlapsedms accumulates the milliseconds spent on the convolutions;
 *   totqrlapsedms accumulates the milliseconds spent on the Householder QR;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs;
 *   totbslapsedms accumulates the milliseconds spent on back substitutions;
 *   totupdlapsedms accumulates the milliseconds spent on updates;
 *   totreslapsedms accumulates the milliseconds spent on residuals;
 *   vrblvl    is the verbose level;
 *   mode      execution mode, 0 (GPU only), 1 (CPU only) or 2 (GPU+CPU).
 *
 * ON RETURN :
 *   tailidx_h is the updated value for tailidx_h;
 *   tailidx_d is the updated value for tailidx_d;
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
 *   resmax    the maximum element of the residual vectors;
 *   zeroQ_h   false if Q was computed on host;
 *   noqr_h    updated flag if ||dx_0|| is zero for the first time on host;
 *   zeroQ_d   false if Q was computed on device;
 *   noqr_d    updated flag if ||dx_0|| is zero for the first time on device;
 *   upidx_h   counts the number of updates skipped by host;
 *   bsidx_h   counts the number of backsubstitutions skipped by host;
 *   upidx_d   counts the number of updates skipped by device;
 *   bsidx_d   counts the number of backsubstitutions skipped by device;
 *   totcnvlapsedms accumulates the milliseconds spent on the convolutions;
 *   totqrlapsedms accumulates the milliseconds spent on the Householder QR;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs;
 *   totbslapsedms accumulates the milliseconds spent on back substitutions;
 *   totupdlapsedms accumulates the milliseconds spent on updates;
 *   totreslapsedms accumulates the milliseconds spent on residuals. */

int dbl_allocate_inoutfunjac
 ( int dim, int deg, int mode,
   double **input_h, double **input_d,
   double ***output_h, double ***output_d,
   double **funval_h, double **funval_d,
   double ***jacval_h, double ***jacval_d );
/*
 * DESCRIPTION :
 *   Allocates work space memory for input, output,
 *   the function values and the value of the Jacobian matrix.
 *
 * ON ENTRY :
 *   dim        dimension of the system;
 *   deg        degree at which the series are truncated;
 *   mode       0 (GPU), 1 (CPU), or 2 (GPU+CPU).
 *
 * ON RETURN :
 *   input_h    input on the host;
 *   input_d    input on the device;
 *   output_h   output on the host;
 *   output_d   output on the device;
 *   funval_h   function values on the host;
 *   funval_d   function values on the device;
 *   jacval_h   values of the Jacobian matrix on the host;
 *   jacval_d   values of the Jacobian matrix on the device. */

int dbl_allocate_rhsqrsol
 ( int dim, int deg, int mode,
   double **rhs_h, double **rhs_d, double **urhs_h, double **urhs_d,
   double **Q_h, double **Q_d, double **R_h, double **R_d,
   double **sol_h, double **sol_d );
/*
 * DESCRIPTION :
 *   Allocates work space memory for the linearized power series system.
 *
 * ON ENTRY :
 *   dim        dimension of the system;
 *   deg        degree at which the series are truncated;
 *   mode       0 (GPU), 1 (CPU), or 2 (GPU+CPU).
 *
 * ON RETURN :
 *   rhs_h      right-hand side on the host;
 *   rhs_d      right-hand side on the device;
 *   urhs_h     updated right-hand side on the host;
 *   urhs_d     updated right-hand side on the device;
 *   Q_h        Q on the host;
 *   Q_d        Q on the device;
 *   R_h        R on the host;
 *   R_d        R on the device;
 *   sol_h      update to the solution on the host;
 *   sol_h      update to the solution of the device. */

void dbl_start_setup
 ( int dim, int deg, double **testsol,
   double **input_h, double **input_d, int mode, int vrblvl );
/*
 * DESCRIPTION :
 *   Given the test solution, defines the start vector.
 *
 * ON ENTRY :
 *   dim        dimension of the solution vectors;
 *   deg        degree at which the series are truncated;
 *   testsol    test solution;
 *   input_h    allocated on host if mode is 1 or 2;
 *   input_d    allocated on device if mode is 0 or 2;
 *   mode       0 (GPU), 1 (CPU), or 2 (GPU+CPU);
 *   vrblvl     is the verbose level.
 *
 * ON RETURN :
 *   input_h    start vector for host if mode is 1 or 2;
 *   input_d    start vector for device if mode is 0 or 2. */

void dbl_column_setup
 ( int dim, int deg, int nbrcol,
   int **nvr, int ***idx, int **rowsA, double ***cff,
   double **testsol, double **mbrhs,
   double **input_h, double **input_d, int mode, int vrblvl );
/*
 * DESCRIPTION :
 *   Defines the test solution and start vectors to run Newton's method
 *   on one or more columns of monomials.
 *
 * ON ENTRY :
 *   dim        dimension of the solution vectors;
 *   deg        degree at which the series are truncated;
 *   nbrcol     number of columns, if 1, then the system is monomial,
 *              otherwise nbrcol columns are expected;
 *   nvr        nvr[i][j] is the number of variables in the j-th monomial
 *              of the i-th column;
 *   idx        idx[i][j] are the indices of the variables in monomial j
 *              of the i-th column;
 *   rowsA      exponents for monomials if only one column;
 *   cff        coefficients of the monomials, if more than one column;
 *   testsol    space for dim pointers;
 *   mbrhs      space for dim pointers;
 *   input_h    allocated on host if mode is 1 or 2;
 *   input_d    allocated on device if mode is 0 or 2;
 *   mode       0 (GPU), 1 (CPU), or 2 (GPU+CPU);
 *   vrblvl     is the verbose level.
 *
 * ON RETURN :
 *   testsol    test solution;
 *   mbrhs      right hand side vector for the test solution;
 *   input_h    start vector for host if mode is 1 or 2;
 *   input_d    start vector for device if mode is 0 or 2. */

void dbl_row_setup
 ( int dim, int deg, int *nbr, int **nvr, int ***idx,
   double **cst, double ***cff, double **testsol, 
   double **input_h, double **input_d,
   double ***output_h, double ***output_d, int mode, int vrblvl );
/*
 * DESCRIPTION :
 *   Defines the test solution and start vectors to run Newton's method
 *   on one or more columns of monomials.
 *
 * ON ENTRY :
 *   dim        dimension of the solution vectors;
 *   deg        degree at which the series are truncated;
 *   nbr        nbr[i] is the number of monomials in the i-th polynomial;
 *   nvr        nvr[i][j] is the number of variables in the j-th monomial
 *              of the i-th polynomial;
 *   idx        idx[i][j] are the indices of the variables in monomial j
 *              of the i-th polynomial;
 *   cst        constants of the polynomials;
 *   cff        coefficients of the monomials;
 *   testsol    space for dim pointers;
 *   input_h    allocated on host if mode is 1 or 2;
 *   input_d    allocated on device if mode is 0 or 2;
 *   output_h   allocated on host if mode is 1 or 2;
 *   output_d   allocated on device if mode is 0 or 2;
 *   mode       0 (GPU), 1 (CPU), or 2 (GPU+CPU);
 *   vrblvl     is the verbose level.
 *
 * ON RETURN :
 *   cst        constant, adjusted for the test solution;
 *   testsol    test solution;
 *   input_h    start vector for host if mode is 1 or 2;
 *   input_d    start vector for device if mode is 0 or 2;
 *   output_h   evaluated test solution if mode is 1;
 *   output_d   evaluated test solution if mode is 0 or 2. */

int dbl_error_testsol
 ( int dim, int deg, int mode,
   double **testsol, double **input_h, double **input_d );
/*
 * DESCRIPTION :
 *   Compares the computed solution against the test solution.
 *
 * ON ENTRY :
 *   dim        dimension of the solution vectors;
 *   deg        degree at which the series are truncated;
 *   mode       0 (GPU), 1 (CPU), or 2 (GPU+CPU).
 *   testsol    series of the test solution;
 *   input_h    results on host if mode is 1 or 2;
 *   input_d    results on device if mode is 0 or 2. */

int test_dbl_column_newton
 ( int szt, int nbt, int dim, int deg, int nbrcol,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac, int **rowsA,
   double dpr, int nbsteps, int mode, int vrblvl );
/*
 * DESCRIPTION :
 *   Runs Newton with real double arithmetic,
 *   on one or more columns of monomials.
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
 *   dpr       damper multiplier for t, should be in (0.0, 1.0];
 *   nbsteps   the number of Newton steps;
 *   mode      the mode of execution, 0 for GPU only, 1 for CPU only,
 *             2 for CPU and GPU;
 *   vrblvl    is the verbose level. */

int test_dbl_row_newton
 ( int szt, int nbt, int dim, int deg, int *nbr, int **nvr, int ***idx,
   double dpr, int nbsteps, int mode, int vrblvl );
/*
 * DESCRIPTION :
 *   Runs Newton with real double arithmetic,
 *   on an indexed polynomial system.
 *
 * ON ENTRY :
 *   szt       size of each tile and block;
 *   nbt       number of tiles and number of blocks;
 *   dim       number of equations and variables in the system;
 *   deg       degree of the power series;
 *   nbr       nbr[i] is the number of monomials in the i-th polynomial;
 *   nvr       nvr[i][j] is the number of variables in the j-th monomial
 *             of the i-th polynomial;
 *   idx       idx[i][j] are the indices of the variables in monomial j
 *             of the i-th polynomial;
 *   dpr       damper multiplier for t, should be in (0.0, 1.0];
 *   nbsteps   the number of Newton steps;
 *   mode      the mode of execution, 0 for GPU only, 1 for CPU only,
 *             2 for CPU and GPU;
 *   vrblvl    is the verbose level. */

#endif
