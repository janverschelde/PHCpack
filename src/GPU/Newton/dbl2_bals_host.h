// The file dbl2_bals.h specifies functions to solve linear systems
// of power series by linearization in double double precision.

#ifndef __dbl2_bals_h__
#define __dbl2_bals_h__

void CPU_dbl2_qrbs_head
 ( int dim, int degp1, double ***mathi, double ***matlo,
   double **rhshi, double **rhslo, double **solhi, double **sollo,
   double **Qhi, double **Qlo, double **Rhi, double **Rlo,
   double *wrkvechi, double *wrkveclo,
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
 *   mathi    degp1 matrices of dimension dim;
 *   matlo    degp1 matrices of dimension dim;
 *   rhshi    degp1 vectors of dimension dim;
 *   rhslo    degp1 vectors of dimension dim;
 *   solhi    space allocated for degp1 vectors of dimension dim;
 *   sollo    space allocated for degp1 vectors of dimension dim;
 *   Qhi      space allocated for a matrix of dimension dim;
 *   Qlo      space allocated for a matrix of dimension dim;
 *   Rhi      space allocated for a matrix of dimension dim;
 *   Rlo      space allocated for a matrix of dimension dim;
 *   wrkvechi is work space allocated for a vector of dimension dim;
 *   wrkveclo is work space allocated for a vector of dimension dim;
 *   zeroQ    if true, then Q is zero and Q must be computed;
 *   noqr     flag if true, then no qr, only when not zeroQ;
 *   vrblvl   the verbose level (0 for silent).
 *
 * ON RETURN :
 *   Qhi      high doubles of the Q of the QR of the Jacobian matrix;
 *   Qlo      low doubles of the Q of the QR of the Jacobian matrix;
 *   Rhi      high doubles of the R in the QR of the Jacobian matrix;
 *   Rlo      low doubles of the R in the QR of the Jacobian matrix;
 *   wrkvechi is work space used to solve the linear system;
 *   wrkveclo is work space used to solve the linear system;
 *   solhi    high doubles of the coefficients of the solution series;
 *   sollo    low doubles of the coefficients of the solution series;
 *   zeroQ    false if Q was computed;
 *   noqr     updated flag if ||dx_0|| is zero for the first time. */

void CPU_cmplx2_qrbs_head
 ( int dim, int degp1,
   double ***matrehi, double ***matrelo, double ***matimhi, double ***matimlo,
   double **rhsrehi, double **rhsrelo, double **rhsimhi, double **rhsimlo,
   double **solrehi, double **solrelo, double **solimhi, double **solimlo,
   double **Qrehi, double **Qrelo, double **Qimhi, double **Qimlo,
   double **Rrehi, double **Rrelo, double **Rimhi, double **Rimlo,
   double *wrkvecrehi, double *wrkvecrelo,
   double *wrkvecimhi, double *wrkvecimlo,
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
 *   matrehi  degp1 matrices of dimension dim;
 *   matrelo  degp1 matrices of dimension dim;
 *   matimhi  degp1 matrices of dimension dim;
 *   matimlo  degp1 matrices of dimension dim;
 *   rhsrehi  degp1 vectors of dimension dim;
 *   rhsrelo  degp1 vectors of dimension dim;
 *   rhsimhi  degp1 vectors of dimension dim;
 *   rhsimlo  degp1 vectors of dimension dim;
 *   solrehi  space allocated for degp1 vectors of dimension dim;
 *   solrelo  space allocated for degp1 vectors of dimension dim;
 *   solimhi  space allocated for degp1 vectors of dimension dim;
 *   solimlo  space allocated for degp1 vectors of dimension dim;
 *   Qrehi    space allocated for a matrix of dimension dim;
 *   Qrelo    space allocated for a matrix of dimension dim;
 *   Qimhi    space allocated for a matrix of dimension dim;
 *   Qimlo    space allocated for a matrix of dimension dim;
 *   Rrehi    space allocated for a matrix of dimension dim;
 *   Rrelo    space allocated for a matrix of dimension dim;
 *   Rimhi    space allocated for a matrix of dimension dim;
 *   Rimlo    space allocated for a matrix of dimension dim;
 *   wrkvecrehi is work space allocated for a vector of dimension dim;
 *   wrkvecrelo is work space allocated for a vector of dimension dim;
 *   wrkvecimhi is work space allocated for a vector of dimension dim;
 *   wrkvecimlo is work space allocated for a vector of dimension dim;
 *   zeroQ    if true, then Q is zero and Q must be computed;
 *   noqr     flag if true, then no qr, only when not zeroQ;
 *   vrblvl   the verbose level (0 for silent).
 *
 * ON RETURN :
 *   Qrehi    high doubles of the real parts of the Q in the QR;
 *   Qrelo    low doubles of the real parts of the Q in the QR;
 *   Qimhi    high doubles of the imaginary parts of the Q in the QR;
 *   Qimlo    low doubles of the imaginary parts of the Q in the QR;
 *   Rrehi    high doubles of the real parts of the R in the QR;
 *   Rrelo    low doubles of the real parts of the R in the QR;
 *   Rimhi    high doubles of the imaginary parts of the R in the QR;
 *   Rimlo    low doubles of the imaginary parts of the R in the QR;
 *   wrkvecrehi is work space used to solve the linear system;
 *   wrkvecrelo is work space used to solve the linear system;
 *   wrkvecimhi is work space used to solve the linear system;
 *   wrkvecimlo is work space used to solve the linear system;
 *   solrehi  high doubles of the real parts of the head of the solution;
 *   solrelo  low doubles of the real parts of the head of the solution;
 *   solimhi  high doubles of the imag parts of the head of the solution;
 *   solimlo  low doubles of the imag parts of the head of the solution;
 *   zeroQ    false if Q was computed;
 *   noqr     updated flag if ||dx_0|| is zero for the first time. */

void CPU_dbl2_qrbs_tail
 ( int dim, int degp1, int tailidx, double ***mathi, double ***matlo,
   double **rhshi, double **rhslo, double **solhi, double **sollo,
   double **Qhi, double **Qlo, double **Rhi, double **Rlo,
   double *wrkvechi, double *wrkveclo, int *upidx, int *bsidx, int *newtail,
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
 *   wrkvechi is work space vector of dimension dim for the substitution;
 *   wrkveclo is work space vector of dimension dim for the substitution;
 *   vrblvl   the verbose level (0 for silent).
 *
 * ON RETURN :
 *   rhshi    updated high doubles of right hand side used as work space;
 *   rhslo    updated low doubles of right hand side used as work space;
 *   solhi    high doubles of all coefficients of the solution series;
 *   sollo    low doubles of all coefficients of the solution series;
 *   upidx    counts the number of updates skipped;
 *   bsidx    counts the number of backsubstitutions skipped;
 *   newtail  the new value for tailidx. */

void CPU_cmplx2_qrbs_tail
 ( int dim, int degp1, int tailidx,
   double ***matrehi, double ***matrelo, double ***matimhi, double ***matimlo,
   double **rhsrehi, double **rhsrelo, double **rhsimhi, double **rhsimlo,
   double **solrehi, double **solrelo, double **solimhi, double **solimlo,
   double **Qrehi, double **Qrelo, double **Qimhi, double **Qimlo,
   double **Rrehi, double **Rrelo, double **Rimhi, double **Rimlo,
   double *wrkvecrehi, double *wrkvecrelo,
   double *wrkvecimhi, double *wrkvecimlo,
   int *upidx, int *bsidx, int *newtail, int vrblvl );
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
 *   matrehi  degp1 matrices of dimension dim;
 *   matrelo  degp1 matrices of dimension dim;
 *   matimhi  degp1 matrices of dimension dim;
 *   matimlo  degp1 matrices of dimension dim;
 *   rhsrehi  degp1 vectors of dimension dim;
 *   rhsrelo  degp1 vectors of dimension dim;
 *   rhsimhi  degp1 vectors of dimension dim;
 *   rhsimlo  degp1 vectors of dimension dim;
 *   solrehi  space allocated for degp1 vectors of dimension dim, with the
 *            high doubles of the real parts of leading coefficients defined;
 *   solrelo  space allocated for degp1 vectors of dimension dim, with the
 *            low doubles of the real parts of leading coefficients defined;
 *   solimhi  space allocated for degp1 vectors of dimension dim, with the
 *            high doubles of the imag parts of leading coefficients defined;
 *   solimlo  space allocated for degp1 vectors of dimension dim, with the
 *            low doubles of the imag parts of leading coefficients defined;
 *   Qrehi    high doubles of the real parts of the Q of the QR;
 *   Qrelo    low doubles of the real parts of the Q of the QR;
 *   Qimhi    high doubles of the imaginary parts of the Q of the QR;
 *   Qimlo    low doubles of the imaginary parts of the Q of the QR;
 *   Rrehi    high doubles of the real parts of the R of the QR;
 *   Rrelo    low doubles of the real parts of the R of the QR;
 *   Rimhi    high doubles of the imaginary parts of the R of the QR;
 *   Rimlo    low doubles of the imaginary parts of the R of the QR;
 *   wrkvecrehi is work space of dimension dim for the substitution;
 *   wrkvecrelo is work space of dimension dim for the substitution;
 *   wrkvecimhi is work space of dimension dim for the substitution;
 *   wrkvecimlo is work space of dimension dim for the substitution;
 *   vrblvl   the verbose level (0 for silent).
 *
 * ON RETURN :
 *   rhsrehi  high doubles of the real parts of updated right hand side;
 *   rhsrelo  low doubles of the real parts of updated right hand side;
 *   rhsimhi  high doubles of the imaginary parts of updated right hand side;
 *   rhsimlo  low doubles of the imaginary parts of updated right hand side;
 *   solrehi  high doubles of the real parts of the solution;
 *   solrelo  low doubles of the real parts of the solution;
 *   solimhi  high doubles of the imaginary parts of the solution;
 *   solimlo  low doubles of the imaginary parts of the solution;
 *   upidx    counts the number of updates skipped;
 *   bsidx    counts the number of backsubstitutions skipped;
 *   newtail  the new value for tailidx. */

void CPU_dbl2_qrbs_solve
 ( int dim, int degp1, int tailidx, double ***mathi, double ***matlo,
   double **rhshi, double **rhslo, double **solhi, double **sollo,
   double **Qhi, double **Qlo, double **Rhi, double **Rlo,
   double *wrkvechi, double *wrkveclo,
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
 *   mathi    degp1 matrices of dimension dim;
 *   matlo    degp1 matrices of dimension dim;
 *   rhshi    degp1 vectors of dimension dim;
 *   rhslo    degp1 vectors of dimension dim;
 *   solhi    space allocated for degp1 vectors of dimension dim;
 *   sollo    space allocated for degp1 vectors of dimension dim;
 *   Qhi      space allocated for a matrix of dimension dim;
 *   Qlo      space allocated for a matrix of dimension dim;
 *   Rhi      space allocated for a matrix of dimension dim;
 *   Rlo      space allocated for a matrix of dimension dim;
 *   wrkvechi has work space allocated for a vector of dimension dim;
 *   wrkveclo has work space allocated for a vector of dimension dim;
 *   zeroQ    if true, then Q is zero and Q must be computed;
 *   noqr     flag if true, then no qr, only when not zeroQ;
 *   vrblvl   the verbose level (0 for silent).
 *
 * ON RETURN :
 *   Qhi      high doubles of the Q of the QR of the Jacobian matrix;
 *   Qlo      low doubles of the Q of the QR of the Jacobian matrix;
 *   Rhi      high doubles of the R in the QR of the Jacobian matrix;
 *   Rlo      low doubles of the R in the QR of the Jacobian matrix;
 *   wrkvechi is work space used to solve the linear system;
 *   wrkveclo is work space used to solve the linear system;
 *   solhi    high double coefficients of the solution series;
 *   sollo    low double coefficients of the solution series;
 *   zeroQ    false if Q was computed;
 *   noqr     updated flag if ||dx_0|| is zero for the first time;
 *   upidx    counts the number of updates skipped;
 *   bsidx    counts the number of backsubstitutions skipped;
 *   newtail  the new value for tailidx. */

void CPU_cmplx2_qrbs_solve
 ( int dim, int degp1, int tailidx,
   double ***matrehi, double ***matrelo, double ***matimhi, double ***matimlo, 
   double **rhsrehi, double **rhsrelo, double **rhsimhi, double **rhsimlo,
   double **solrehi, double **solrelo, double **solimhi, double **solimlo,
   double **Qrehi, double **Qrelo, double **Qimhi, double **Qimlo,
   double **Rrehi, double **Rrelo, double **Rimhi, double **Rimlo,
   double *wrkvecrehi, double *wrkvecrelo,
   double *wrkvecimhi, double *wrkvecimlo,
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
 *   matrehi  degp1 matrices of dimension dim;
 *   matrelo  degp1 matrices of dimension dim;
 *   matimhi  degp1 matrices of dimension dim;
 *   matimlo  degp1 matrices of dimension dim;
 *   rhsrehi  degp1 vectors of dimension dim;
 *   rhsrelo  degp1 vectors of dimension dim;
 *   rhsimhi  degp1 vectors of dimension dim;
 *   rhsimlo  degp1 vectors of dimension dim;
 *   solrehi  space allocated for degp1 vectors of dimension dim;
 *   solrelo  space allocated for degp1 vectors of dimension dim;
 *   solimhi  space allocated for degp1 vectors of dimension dim;
 *   solimlo  space allocated for degp1 vectors of dimension dim;
 *   Qrehi    space allocated for a matrix of dimension dim;
 *   Qrelo    space allocated for a matrix of dimension dim;
 *   Qimhi    space allocated for a matrix of dimension dim;
 *   Qimlo    space allocated for a matrix of dimension dim;
 *   Rrehi    space allocated for a matrix of dimension dim;
 *   Rrelo    space allocated for a matrix of dimension dim;
 *   Rimhi    space allocated for a matrix of dimension dim;
 *   Rimlo    space allocated for a matrix of dimension dim;
 *   wrkvecrehi has work space allocated for a vector of dimension dim;
 *   wrkvecrelo has work space allocated for a vector of dimension dim;
 *   wrkvecimhi has work space allocated for a vector of dimension dim;
 *   wrkvecimlo has work space allocated for a vector of dimension dim;
 *   zeroQ    if true, then Q is zero and Q must be computed;
 *   noqr     flag if true, then no qr, only when not zeroQ;
 *   vrblvl   the verbose level (0 for silent).
 *
 * ON RETURN :
 *   Qrehi    high doubles of the real parts of the Q in the QR;
 *   Qrelo    low doubles of the real parts of the Q in the QR;
 *   Qimhi    high doubles of the imaginary parts of the Q in the QR;
 *   Qimlo    low doubles of the imaginary parts of the Q in the QR;
 *   Rrehi    high doubles of the real parts of the R in the QR;
 *   Rrelo    low doubles of the real parts of the R in the QR;
 *   Rimhi    high doubles of the imaginary parts of the R in the QR;
 *   Rimlo    low doubles of the imaginary parts of the R in the QR;
 *   wrkvecrehi is work space used to solve the linear systems;
 *   wrkvecrelo is work space used to solve the linear systems;
 *   wrkvecimhi is work space used to solve the linear systems;
 *   wrkvecimlo is work space used to solve the linear systems;
 *   solrehi  high doubles of the real parts of the solution;
 *   solrelo  low doubles of the real parts of the solution;
 *   solimhi  high doubles of the imaginary parts of the solution;
 *   solimlo  low doubles of the imaginary parts of the solution;
 *   zeroQ    false if Q was computed;
 *   noqr     updated flag if ||dx_0|| is zero for the first time;
 *   upidx    counts the number of updates skipped;
 *   bsidx    counts the number of backsubstitutions skipped;
 *   newtail  the new value for tailidx. */

void CPU_dbl2_linear_residue
 ( int dim, int degp1, int tailidx, double ***mathi, double ***matlo,
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
 *   tailidx  the index of the start of the update in the tail;
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

void CPU_cmplx2_linear_residue
 ( int dim, int degp1, int tailidx,
   double ***matrehi, double ***matrelo, double ***matimhi, double ***matimlo,
   double **rhsrehi, double **rhsrelo, double **rhsimhi, double **rhsimlo, 
   double **solrehi, double **solrelo, double **solimhi, double **solimlo,
   double **resvecrehi, double **resvecrelo,
   double **resvecimhi, double **resvecimlo,
   double *resmaxhi, double *resmaxlo, int vrblvl );
/*
 * DESCRIPTION :
 *   Computes the residual of the linear power series system.
 *
 * ON ENTRY :
 *   dim      the dimension of the matrices and vectors;
 *   degp1    degree plus one, the size of the matrix system;
 *   tailidx  the index of the start of the update in the tail;
 *   matrehi  degp1 matrices of dimension dim;
 *   matrelo  degp1 matrices of dimension dim;
 *   matimhi  degp1 matrices of dimension dim;
 *   matimlo  degp1 matrices of dimension dim;
 *   rhsrehi  degp1 right hand side vectors of dimension dim;
 *   rhsrelo  degp1 right hand side vectors of dimension dim;
 *   rhsimhi  degp1 right hand side vectors of dimension dim;
 *   rhsimlo  degp1 right hand side vectors of dimension dim;
 *   solrehi  degp1 solution vectors of dimension dim;
 *   solrelo  degp1 solution vectors of dimension dim;
 *   solimhi  degp1 solution vectors of dimension dim;
 *   solimlo  degp1 solution vectors of dimension dim;
 *   resvecrehi has space for the residual power series;
 *   resvecrelo has space for the residual power series;
 *   resvecimhi has space for the residual power series;
 *   resvecimlo has space for the residual power series;
 *   vrblvl   is the verbose level.
 *
 * ON RETURN :
 *   resvecrehi are the high doubles of the real parts of the residual series;
 *   resvecrelo are the low doubles of the real parts of the residual series;
 *   resvecimhi are the high doubles of the imag parts the residual series;
 *   resvecimlo are the low doubles of the imag parts the residual series;
 *   resmaxhi is the high double of the max norm of the residual series;
 *   resmaxlo is the low double of the max norm of the residual series. */

#endif
