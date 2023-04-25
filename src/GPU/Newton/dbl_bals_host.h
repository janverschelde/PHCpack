// The file dbl_bals_host.h specifies functions to solve linear systems
// of power series by linearization in double precision.

#ifndef __dbl_bals_host_h__
#define __dbl_bals_host_h__

void CPU_dbl_qrbs_head
 ( int dim, int degp1, double ***mat, double **rhs, double **sol,
   double **Q, double **R, double *wrkvec,
   bool *zeroQ, bool *noqr, int vrblvl );
/*
 * DESCRIPTION :
 *   Computes the leading terms of the power series solution
 *   to a linear system of power series, in linearized format,
 *   using back substitution after a QR factorization, on real data.
 *
 * ON ENTRY :
 *   dim      the dimension of the matrices and vectors;
 *   degp1    degree plus one, the size of the matrix system;
 *   mat      degp1 matrices of dimension dim;
 *   rhs      degp1 vectors of dimension dim;
 *   sol      space allocated for degp1 vectors of dimension dim;
 *   Q        space allocated for a matrix of dimension dim;
 *   R        space allocated for a matrix of dimension dim;
 *   wrkvec   work space allocated for a vector of dimension dim;
 *   zeroQ    if true, then Q is zero and Q must be computed;
 *   noqr     flag if true, then no qr, only when not zeroQ;
 *   vrblvl   the verbose level (0 for silent).
 *
 * ON RETURN :
 *   Q        the Q in a QR factorization of the Jacobian matrix;
 *   R        the R in a QR factorization of the Jacobian matrix;
 *   wrkvec   work space used to solve the linear system;
 *   sol      the coefficients of the solution series;
 *   zeroQ    false if Q was computed;
 *   noqr     updated flag if ||dx_0|| is zero for the first time. */

void CPU_cmplx_qrbs_head
 ( int dim, int degp1, double ***matre, double ***matim,
   double **rhsre, double **rhsim, double **solre, double **solim,
   double **Qre, double **Qim, double **Rre, double **Rim,
   double *wrkvecre, double *wrkvecim,
   bool *zeroQ, bool *noqr, int vrblvl );
/*
 * DESCRIPTION :
 *   Computes the leading terms of the power series solution
 *   to a linear system of power series, in linearized format,
 *   using back substitution after a QR factorization, on complex data.
 *
 * ON ENTRY :
 *   dim      the dimension of the matrices and vectors;
 *   degp1    degree plus one, the size of the matrix system;
 *   matre    degp1 matrices of dimension dim;
 *   matim    degp1 matrices of dimension dim;
 *   rhsre    degp1 vectors of dimension dim;
 *   rhsim    degp1 vectors of dimension dim;
 *   solre    space allocated for degp1 vectors of dimension dim;
 *   solim    space allocated for degp1 vectors of dimension dim;
 *   Qre      space allocated for a matrix of dimension dim;
 *   Qim      space allocated for a matrix of dimension dim;
 *   Rre      space allocated for a matrix of dimension dim;
 *   Rim      space allocated for a matrix of dimension dim;
 *   wrkvecre is work space allocated for a vector of dimension dim;
 *   wrkvecim is work space allocated for a vector of dimension dim;
 *   zeroQ    if true, then Q is zero and Q must be computed;
 *   noqr     flag if true, then no qr, only when not zeroQ;
 *   vrblvl   the verbose level (0 for silent).
 *
 * ON RETURN :
 *   Qre      real parts of the Q in the QR of the Jacobian;
 *   Qim      imaginary parts of the Q in the QR of the Jacobian;
 *   Rre      real parts of the R in the QR of the Jacobian;
 *   Rim      imaginary parts of the R in the QR of the Jacobian;
 *   wrkvecre is work space used to solve the linear system;
 *   wrkvecim is work space used to solve the linear system;
 *   solre    real parts of the head term of the solution series;
 *   solim    imaginary parts of the head term of the solution series;
 *   zeroQ    false if Q was computed;
 *   noqr     updated flag if ||dx_0|| is zero for the first time. */

void CPU_dbl_qrbs_tail
 ( int dim, int degp1, int tailidx, double ***mat, double **rhs, double **sol,
   double **Q, double **R, double *wrkvec, int *upidx, int *bsidx,
   int *newtail, int vrblvl );
/*
 * DESCRIPTION :
 *   Computes the trailing terms of the power series solution
 *   to a linear system of power series, in linearized format,
 *   applying substitution given a QR factorization.
 *
 * ON ENTRY :
 *   dim      the dimension of the matrices and vectors;
 *   degp1    degree plus one, the size of the matrix system;
 *   tailidx  the index of the start of the update in the tail;
 *   mat      degp1 matrices of dimension dim;
 *   rhs      degp1 vectors of dimension dim;
 *   sol      space allocated for degp1 vectors of dimension dim,
 *            with the leading coefficients defined;
 *   Q        has the Q of the QR factorization;
 *   R        has the R of the QR factorization;
 *   wrkvec   work space vector of dimension dim for the substitution;
 *   vrblvl   the verbose level (0 for silent).
 *
 * ON RETURN :
 *   rhs      updated right hand side used as work space;
 *   sol      all coefficients of the solution series;
 *   upidx    counts the number of updates skipped;
 *   bsidx    counts the number of backsubstitutions skipped;
 *   newtail  the new value for tailidx. */

void CPU_cmplx_qrbs_tail
 ( int dim, int degp1, int tailidx, double ***matre, double ***matim,
   double **rhsre, double **rhsim, double **solre, double **solim,
   double **Qre, double **Qim, double **Rre, double **Rim,
   double *wrkvecre, double *wrkvecim, int *upidx, int *bsidx, int *newtail,
   int vrblvl );
/*
 * DESCRIPTION :
 *   Computes the trailing terms of the power series solution
 *   to a linear system of power series, in linearized format,
 *   applying substitution given a QR factorization.
 *
 * ON ENTRY :
 *   dim      the dimension of the matrices and vectors;
 *   degp1    degree plus one, the size of the matrix system;
 *   tailidx  the index of the start of the update in the tail;
 *   matre    degp1 matrices of dimension dim;
 *   matim    degp1 matrices of dimension dim;
 *   rhsre    degp1 vectors of dimension dim;
 *   rhsim    degp1 vectors of dimension dim;
 *   solre    space allocated for degp1 vectors of dimension dim,
 *            with the real parts of leading coefficients defined;
 *   solim    space allocated for degp1 vectors of dimension dim,
 *            with the imaginary parts of leading coefficients defined;
 *   Qre      real parts of the Q of the QR factorization;
 *   Qim      imaginary parts of the Q of the QR factorization;
 *   Rre      real parts of the R of the QR factorization;
 *   Rim      imaginary parts of the R of the QR factorization;
 *   wrkvecre is work space of dimension dim for the substitution;
 *   wrkvecim is work space of dimension dim for the substitution;
 *   vrblvl   the verbose level (0 for silent).
 *
 * ON RETURN :
 *   rhsre    real parts of updated right hand side;
 *   rhsim    imaginary parts of updated right hand side;
 *   solre    real parts of all coefficients of the solution;
 *   solim    imaginary parts of all coefficients of the solution;
 *   upidx    counts the number of updates skipped;
 *   bsidx    counts the number of backsubstitutions skipped;
 *   newtail  the new value for tailidx. */

void CPU_dbl_qrbs_solve
 ( int dim, int degp1, int tailidx, double ***mat, double **rhs, double **sol,
   double **Q, double **R, double *wrkvec,
   bool *zeroQ, bool *noqr, int *upidx, int *bsidx, int *newtail,
   int vrblvl );
/*
 * DESCRIPTION :
 *   Solves a linear system of power series, in linearized format,
 *   using QR factorization and substitutions.
 *
 * ON ENTRY :
 *   dim      the dimension of the matrices and vectors;
 *   degp1    degree plus one, the size of the matrix system;
 *   tailidx  the index of the start of the update in the tail;
 *   mat      degp1 matrices of dimension dim;
 *   rhs      degp1 vectors of dimension dim;
 *   sol      space allocated for degp1 vectors of dimension dim;
 *   Q        space allocated for a matrix of dimension dim;
 *   R        space allocated for a matrix of dimension dim;
 *   wrkvec   work space allocated for a vector of dimension dim;
 *   zeroQ    if true, then Q is zero and Q must be computed;
 *   noqr     flag if true, then no qr, only when not zeroQ;
 *   vrblvl   the verbose level (0 for silent).
 *
 * ON RETURN :
 *   Q        the Q in a QR factorization of the Jacobian matrix;
 *   R        the R in a QR factorization of the Jacobian matrix;
 *   wrkvec   work space used to solve the linear system;
 *   sol      the coefficients of the solution series;
 *   zeroQ    false if Q was computed;
 *   noqr     updated flag if ||dx_0|| is zero for the first time;
 *   upidx    counts the number of updates skipped;
 *   bsidx    counts the number of backsubstitutions skipped;
 *   newtail  the new value for tailidx. */

void CPU_cmplx_qrbs_solve
 ( int dim, int degp1, int tailidx, double ***matre, double ***matim, 
   double **rhsre, double **rhsim, double **solre, double **solim,
   double **Qre, double **Qim, double **Rre, double **Rim,
   double *wrkvecre, double *wrkvecim,
   bool *zeroQ, bool *noqr, int *upidx, int *bsidx, int *newtail,
   int vrblvl );
/*
 * DESCRIPTION :
 *   Solves a linear system of power series, in linearized format,
 *   using QR factorization and substitutions.
 *
 * ON ENTRY :
 *   dim      the dimension of the matrices and vectors;
 *   degp1    degree plus one, the size of the matrix system;
 *   tailidx  the index of the start of the update in the tail;
 *   matre    degp1 matrices of dimension dim;
 *   matim    degp1 matrices of dimension dim;
 *   rhsre    degp1 vectors of dimension dim;
 *   rhsim    degp1 vectors of dimension dim;
 *   solre    space allocated for degp1 vectors of dimension dim;
 *   solim    space allocated for degp1 vectors of dimension dim;
 *   Qre      space allocated for a matrix of dimension dim;
 *   Qim      space allocated for a matrix of dimension dim;
 *   Rre      space allocated for a matrix of dimension dim;
 *   Rim      space allocated for a matrix of dimension dim;
 *   wrkvecre has work space allocated for a vector of dimension dim;
 *   wrkvecim has work space allocated for a vector of dimension dim;
 *   zeroQ    if true, then Q is zero and Q must be computed;
 *   noqr     flag if true, then no qr, only when not zeroQ;
 *   vrblvl   the verbose level (0 for silent).
 *
 * ON RETURN :
 *   Qre      real parts of the Q in the QR of the Jacobian;
 *   Qim      imaginary parts of the Q in the QR of the Jacobian;
 *   Rre      real parts of the R in the QR of the Jacobian;
 *   Rim      imaginary parts of the R in the QR of the Jacobian;
 *   wrkvecre is work space used to solve the linear systems;
 *   wrkvecim is work space used to solve the linear systems;
 *   solre    real parts of the coefficients of the solution;
 *   solim    imaginary parts of the coefficients of the solution;
 *   zeroQ    false if Q was computed;
 *   noqr     updated flag if ||dx_0|| is zero for the first time;
 *   upidx    counts the number of updates skipped;
 *   bsidx    counts the number of backsubstitutions skipped;
 *   newtail  the new value for tailidx. */

void CPU_dbl_linear_residue
 ( int dim, int degp1, int tailidx, double ***mat, double **rhs, double **sol,
   double **resvec, double *resmax, int vrblvl );
/*
 * DESCRIPTION :
 *   Computes the residual of the linear power series system.
 *
 * ON ENTRY :
 *   dim      the dimension of the matrices and vectors;
 *   degp1    degree plus one, the size of the matrix system;
 *   tailidx  the index of the start of the update in the tail;
 *   mat      degp1 matrices of dimension dim;
 *   rhs      degp1 right hand side vectors of dimension dim;
 *   sol      degp1 solution vectors of dimension dim;
 *   resvec   space for the residual power series;
 *   vrblvl   is the verbose level.
 *
 * ON RETURN :
 *   resvec   the residual power series;
 *   resmax   maximum component of the residual power series. */

void CPU_cmplx_linear_residue
 ( int dim, int degp1, int tailidx, double ***matre, double ***matim,
   double **rhsre, double **rhsim, double **solre, double **solim,
   double **resvecre, double **resvecim, double *resmax, int vrblvl );
/*
 * DESCRIPTION :
 *   Computes the residual of the linear power series system.
 *
 * ON ENTRY :
 *   dim      the dimension of the matrices and vectors;
 *   degp1    degree plus one, the size of the matrix system;
 *   tailidx  the index of the start of the update in the tail;
 *   matre    degp1 matrices of dimension dim;
 *   matim    degp1 matrices of dimension dim;
 *   rhsre    degp1 right hand side vectors of dimension dim;
 *   rhsim    degp1 right hand side vectors of dimension dim;
 *   solre    degp1 solution vectors of dimension dim;
 *   solim    degp1 solution vectors of dimension dim;
 *   resvecre has space for the residual power series;
 *   resvecim has space for the residual power series;
 *   vrblvl   is the verbose level.
 *
 * ON RETURN :
 *   resvecre are the real parts of the residual power series;
 *   resvecim are the imaginary parts the residual power series;
 *   resmax   max norm of the residual power series. */

#endif
