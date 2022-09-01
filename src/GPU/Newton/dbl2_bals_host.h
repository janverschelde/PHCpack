// The file dbl2_bals.h specifies functions to solve linear systems
// of power series by linearization in double double precision.

#ifndef __dbl2_bals_h__
#define __dbl2_bals_h__

void CPU_dbl2_lusb_head
 ( int dim, int degp1, double ***mathi, double ***matlo, 
   double **rhshi, double **rhslo, double **solhi, double **sollo,
   double **wrkmathi, double **wrkmatlo, double *wrkvechi, double *wrkveclo,
   int *pivots, int vrblvl );
/*
 * DESCRIPTION :
 *   Computes the leading terms of the power series solution
 *   to a linear system of power series, in linearized format,
 *   using substitution after a LU factorization.
 *
 * ON ENTRY :
 *   dim      the dimension of the matrices and vectors;
 *   degp1    degree plus one, the size of the matrix system;
 *   mathi    degp1 matrices of dimension dim;
 *   matlo    degp1 matrices of dimension dim;
 *   rhshi    degp1 vectors of dimension dim;
 *   rhslo    degp1 vectors of dimension dim;
 *   solhi    space allocated for degp1 vectors of dimension dim;
 *   sollo    space allocated for degp1 vectors of dimension dim;
 *   wrkmathi is work space allocated for a matrix of dimension dim;
 *   wrkmatlo is work space allocated for a matrix of dimension dim;
 *   wrkvechi is work space allocated for a vector of dimension dim;
 *   wrkveclo is work space allocated for a vector of dimension dim;
 *   pivots   space for dim integers to store the pivots;
 *   vrblvl   the verbose level (0 for silent).
 *
 * ON RETURN :
 *   wrkmathi are high doubles of the LU factorization of the Jacobian matrix;
 *   wrkmatlo are low doubles of the LU factorization of the Jacobian matrix;
 *   pivots   pivots used in the LU factorization of the Jacobian;
 *   solhi    the high doubles of the solution series;
 *   sollo    the low doubles of the solution series. */

void CPU_dbl2_qrbs_head
 ( int dim, int degp1, double ***mathi, double ***matlo,
   double **rhshi, double **rhslo, double **solhi, double **sollo,
   double **wrkmathi, double **wrkmatlo, double **Qhi, double **Qlo,
   double **Rhi, double **Rlo, double *wrkvechi, double *wrkveclo,
   int vrblvl );
/*
 * DESCRIPTION :
 *   Computes the leading terms of the power series solution
 *   to a linear system of power series, in linearized format,
 *   using back substitution after a QR factorization.
 *
 * ON ENTRY :
 *   dim      the dimension of the matrices and vectors;
 *   degp1    degree plus one, the size of the matrix system;
 *   mathi    degp1 matrices of dimension dim;
 *   matlo    degp1 matrices of dimension dim;
 *   rhshi    degp1 vectors of dimension dim;
 *   rhslo    degp1 vectors of dimension dim;
 *   solhi    space allocated for degp1 vectors of dimension dim;
 *   sollo    space allocated for degp1 vectors of dimension dim;
 *   wrkmathi is work space allocated for a matrix of dimension dim;
 *   wrkmatlo is work space allocated for a matrix of dimension dim;
 *   wrkvechi is work space allocated for a vector of dimension dim;
 *   wrkveclo is work space allocated for a vector of dimension dim;
 *   Qhi      space allocated for a matrix of dimension dim;
 *   Qlo      space allocated for a matrix of dimension dim;
 *   Rhi      space allocated for a matrix of dimension dim;
 *   Rlo      space allocated for a matrix of dimension dim;
 *   wrkvechi is work space allocated for a vector of dimension dim;
 *   wrkveclo is work space allocated for a vector of dimension dim;
 *   vrblvl   the verbose level (0 for silent).
 *
 * ON RETURN :
 *   wrkmathi contains a copy of the Jacobian matrix in mathi[0];
 *   wrkmatlo contains a copy of the Jacobian matrix in matlo[0];
 *   Qhi      high doubles of the Q of the QR of the Jacobian matrix;
 *   Qlo      low doubles of the Q of the QR of the Jacobian matrix;
 *   Rhi      high doubles of the R in the QR of the Jacobian matrix;
 *   Rhi      low doubles of the R in the QR of the Jacobian matrix;
 *   wrkvechi is work space used to solve the linear system;
 *   wrkveclo is work space used to solve the linear system;
 *   solhi    high doubles of the coefficients of the solution series;
 *   sollo    low doubles of the coefficients of the solution series. */

void CPU_dbl2_lusb_tail
 ( int dim, int degp1, double ***mathi, double ***matlo,
   double **rhshi, double **rhslo, double **solhi, double **sollo,
   double **wrkmathi, double **wrkmatlo, int *pivots, int vrblvl );
/*
 * DESCRIPTION :
 *   Computes the trailing terms of the power series solution
 *   to a linear system of power series, in linearized format,
 *   applying substitution given a LU factorization.
 *
 * ON ENTRY :
 *   dim      the dimension of the matrices and vectors;
 *   degp1    degree plus one, the size of the matrix system;
 *   mathi    degp1 matrices of dimension dim;
 *   matlo    degp1 matrices of dimension dim;
 *   rhshi    degp1 vectors of dimension dim;
 *   rhslo    degp1 vectors of dimension dim;
 *   solhi    space allocated for degp1 vectors of dimension dim,
 *            with the high doubles of the leading coefficients defined;
 *   sollo    space allocated for degp1 vectors of dimension dim,
 *            with the low doubles of the leading coefficients defined;
 *   wrkmathi has the high part of the LU factorization of the Jacobian matrix;
 *   wrkmatlo has the low part of the LU factorization of the Jacobian matrix;
 *   pivots   dim integers store the pivots of the LU factorization;
 *   vrblvl   the verbose level (0 for silent).
 *
 * ON RETURN :
 *   rhshi    updated right hand side used as work space;
 *   rhslo    updated right hand side used as work space;
 *   solhi    all high double coefficients of the solution series;
 *   sollo    all low double coefficients of the solution series. */

void CPU_dbl2_qrbs_tail
 ( int dim, int degp1, double ***mathi, double ***matlo,
   double **rhshi, double **rhslo, double **solhi, double **sollo,
   double **Qhi, double **Qlo, double **Rhi, double **Rlo,
   double *wrkvechi, double *wrkveclo, int vrblvl );
/*
 * DESCRIPTION :
 *   Computes the trailing terms of the power series solution
 *   to a linear system of power series, in linearized format,
 *   applying substitution given a QR factorization.
 *
 * ON ENTRY :
 *   dim      the dimension of the matrices and vectors;
 *   degp1    degree plus one, the size of the matrix system;
 *   mathi    degp1 matrices of dimension dim;
 *   matlo    degp1 matrices of dimension dim;
 *   rhshi    degp1 vectors of dimension dim;
 *   rhslo    degp1 vectors of dimension dim;
 *   solhi    space allocated for degp1 vectors of dimension dim,
 *            with the leading high doubles defined;
 *   sollo    space allocated for degp1 vectors of dimension dim,
 *            with the leading low doubles defined;
 *   Qhi      high doubles of the Q of the QR factorization;
 *   Qlo      low doubles of the Q of the QR factorization;
 *   Rhi      high doubles of the R of the QR factorization;
 *   Rlo      low doubles of the R of the QR factorization;
 *   wrkvec   work space vector of dimension dim for the substitution;
 *   wrkvec   work space vector of dimension dim for the substitution;
 *   vrblvl   the verbose level (0 for silent).
 *
 * ON RETURN :
 *   rhshi    updated high doubles of right hand side used as work space;
 *   rhslo    updated low doubles of right hand side used as work space;
 *   solhi    high doubles of all coefficients of the solution series;
 *   sollo    low doubles of all coefficients of the solution series. */

void CPU_dbl2_lusb_solve
 ( int dim, int degp1, double ***mathi, double ***matlo,
   double **rhshi, double **rhslo, double **solhi, double **sollo,
   double **wrkmathi, double **wrkmatlo, double *wrkvechi, double *wrkveclo,
   int *pivots, int vrblvl );
/*
 * DESCRIPTION :
 *   Solves a linear system of power series, in linearized format,
 *   using LU factorization and substitutions.
 *
 * ON ENTRY :
 *   dim      the dimension of the matrices and vectors;
 *   degp1    degree plus one, the size of the matrix system;
 *   mathi    degp1 matrices of dimension dim;
 *   matlo    degp1 matrices of dimension dim;
 *   rhshi    degp1 vectors of dimension dim;
 *   rhslo    degp1 vectors of dimension dim;
 *   solhi    space allocated for degp1 vectors of dimension dim;
 *   sollo    space allocated for degp1 vectors of dimension dim;
 *   wrkmathi has work space allocated for a matrix of dimension dim;
 *   wrkmatlo has work space allocated for a matrix of dimension dim;
 *   wrkvechi has work space allocated for a vector of dimension dim;
 *   wrkveclo has work space allocated for a vector of dimension dim;
 *   pivots   space for dim integers to store the pivots;
 *   vrblvl   the verbose level (0 for silent).
 *
 * ON RETURN :
 *   wrkmathi has the high part of the LU factorization of the Jacobian matrix;
 *   wrkmatlo has the low part of the LU factorization of the Jacobian matrix;
 *   pivots   pivots used in the LU factorization of the Jacobian;
 *   solhi    the high double coefficients of the solution series;
 *   sollo    the low double coefficients of the solution series. */

void CPU_dbl2_qrbs_solve
 ( int dim, int degp1, double ***mathi, double ***matlo,
   double **rhshi, double **rhslo, double **solhi, double **sollo,
   double **wrkmathi, double **wrkmatlo,
   double **Qhi, double **Qlo, double **Rhi, double **Rlo,
   double *wrkvechi, double *wrkveclo, int vrblvl );
/*
 * DESCRIPTION :
 *   Solves a linear system of power series, in linearized format,
 *   using QR factorization and substitutions.
 *
 * ON ENTRY :
 *   dim      the dimension of the matrices and vectors;
 *   degp1    degree plus one, the size of the matrix system;
 *   mathi    degp1 matrices of dimension dim;
 *   matlo    degp1 matrices of dimension dim;
 *   rhshi    degp1 vectors of dimension dim;
 *   rhslo    degp1 vectors of dimension dim;
 *   solhi    space allocated for degp1 vectors of dimension dim;
 *   sollo    space allocated for degp1 vectors of dimension dim;
 *   wrkmathi has work space allocated for a matrix of dimension dim;
 *   wrkmatlo has work space allocated for a matrix of dimension dim;
 *   Qhi      space allocated for a matrix of dimension dim;
 *   Qlo      space allocated for a matrix of dimension dim;
 *   Rhi      space allocated for a matrix of dimension dim;
 *   Rlo      space allocated for a matrix of dimension dim;
 *   wrkvechi has work space allocated for a vector of dimension dim;
 *   wrkveclo has work space allocated for a vector of dimension dim;
 *   vrblvl   the verbose level (0 for silent).
 *
 * ON RETURN :
 *   wrkmathi contains a copy of the Jacobian matrix in mathi[0];
 *   wrkmatlo contains a copy of the Jacobian matrix in matlo[0];
 *   Qhi      high doubles of the Q of the QR of the Jacobian matrix;
 *   Qlo      low doubles of the Q of the QR of the Jacobian matrix;
 *   Rhi      high doubles of the R in the QR of the Jacobian matrix;
 *   Rhi      low doubles of the R in the QR of the Jacobian matrix;
 *   wrkvechi is work space used to solve the linear system;
 *   wrkveclo is work space used to solve the linear system;
 *   solhi    high double coefficients of the solution series;
 *   sollo    low double coefficients of the solution series. */

void CPU_dbl2_linear_residue
 ( int dim, int degp1, double ***mathi, double ***matlo,
   double **rhshi, double **rhslo, double **solhi, double **sollo,
   double **resvechi, double **resveclo, double *resmaxhi, double *resmaxlo,
   int vrblvl );
/*
 * DESCRIPTION :
 *   Computes the residual of the linear power series system.
 *
 * ON ENTRY :
 *   dim      the dimension of the matrices and vectors;
 *   degp1    degree plus one, the size of the matrix system;
 *   mathi    degp1 matrices of dimension dim;
 *   matlo    degp1 matrices of dimension dim;
 *   rhshi    degp1 right hand side vectors of dimension dim;
 *   rhslo    degp1 right hand side vectors of dimension dim;
 *   solhi    degp1 solution vectors of dimension dim;
 *   sollo    degp1 solution vectors of dimension dim;
 *   resvechi has space for the residual power series;
 *   resveclo has space for the residual power series;
 *   vrblvl   is the verbose level.
 *
 * ON RETURN :
 *   resvechi are the high doubles of the residual power series;
 *   resveclo are the low doubles of the residual power series;
 *   resmaxhi is the high double of the maximum component in resvec;
 *   resmaxlo is the low double of the maximum component in resvec. */

#endif
