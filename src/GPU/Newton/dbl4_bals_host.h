// The file dbl4_bals_host.h specifies functions to solve linear systems
// of power series by linearization in quad double precision.

#ifndef __dbl4_bals_host_h__
#define __dbl4_bals_host_h__

void CPU_dbl4_qrbs_head
 ( int dim, int degp1,
   double ***mathihi, double ***matlohi, double ***mathilo, double ***matlolo,
   double **rhshihi, double **rhslohi, double **rhshilo, double **rhslolo,
   double **solhihi, double **sollohi, double **solhilo, double **sollolo,
   double **Qhihi, double **Qlohi, double **Qhilo, double **Qlolo,
   double **Rhihi, double **Rlohi, double **Rhilo, double **Rlolo, 
   double *wrkvechihi, double *wrkveclohi,
   double *wrkvechilo, double *wrkveclolo,
   bool *zeroQ, bool *noqr, int vrblvl );
/*
 * DESCRIPTION :
 *   Computes the leading terms of the power series solution
 *   to a linear system of power series, in linearized format,
 *   using back substitution after a QR factorization.
 *
 * ON ENTRY :
 *   dim      the dimension of the matrices and vectors;
 *   degp1    degree plus one, the size of the matrix system;
 *   mathihi  degp1 matrices of dimension dim;
 *   matlohi  degp1 matrices of dimension dim;
 *   mathilo  degp1 matrices of dimension dim;
 *   matlolo  degp1 matrices of dimension dim;
 *   rhshihi  degp1 vectors of dimension dim;
 *   rhslohi  degp1 vectors of dimension dim;
 *   rhshilo  degp1 vectors of dimension dim;
 *   rhslolo  degp1 vectors of dimension dim;
 *   solhihi  space allocated for degp1 vectors of dimension dim;
 *   sollohi  space allocated for degp1 vectors of dimension dim;
 *   solhilo  space allocated for degp1 vectors of dimension dim;
 *   sollolo  space allocated for degp1 vectors of dimension dim;
 *   wrkvechihi is work space allocated for a vector of dimension dim;
 *   wrkveclohi is work space allocated for a vector of dimension dim;
 *   wrkvechilo is work space allocated for a vector of dimension dim;
 *   wrkveclolo is work space allocated for a vector of dimension dim;
 *   Qhihi    space allocated for a matrix of dimension dim;
 *   Qlohi    space allocated for a matrix of dimension dim;
 *   Qhilo    space allocated for a matrix of dimension dim;
 *   Qlolo    space allocated for a matrix of dimension dim;
 *   Rhihi    space allocated for a matrix of dimension dim;
 *   Rlohi    space allocated for a matrix of dimension dim;
 *   Rhilo    space allocated for a matrix of dimension dim;
 *   Rlolo    space allocated for a matrix of dimension dim;
 *   wrkvechihi is work space allocated for a vector of dimension dim;
 *   wrkveclohi is work space allocated for a vector of dimension dim;
 *   wrkvechilo is work space allocated for a vector of dimension dim;
 *   wrkveclolo is work space allocated for a vector of dimension dim;
 *   zeroQ    if true, then Q is zero and Q must be computed;
 *   noqr     flag if true, then no qr, only when not zeroQ;
 *   vrblvl   the verbose level (0 for silent).
 *
 * ON RETURN :
 *   Qhihi    highest doubles of the Q of the QR of the Jacobian matrix;
 *   Qlohi    second highest doubles of the Q of the QR of the Jacobian matrix;
 *   Qhilo    second lowest doubles of the Q of the QR of the Jacobian matrix;
 *   Qlolo    lowest doubles of the Q of the QR of the Jacobian matrix;
 *   Rhihi    highest doubles of the R in the QR of the Jacobian matrix;
 *   Rlohi    second highest doubles of the R in the QR of the Jacobian matrix;
 *   Rhilo    second lowest doubles of the R in the QR of the Jacobian matrix;
 *   Rlolo    lowest doubles of the R in the QR of the Jacobian matrix;
 *   wrkvechihi is work space used to solve the linear system;
 *   wrkveclohi is work space used to solve the linear system;
 *   wrkvechilo is work space used to solve the linear system;
 *   wrkveclolo is work space used to solve the linear system;
 *   solhihi  highest doubles of the coefficients of the solution;
 *   sollohi  second highest doubles of the coefficients of the solution;
 *   solhilo  second lowest doubles of the coefficients of the solution;
 *   sollolo  lowest doubles of the coefficients of the solution;
 *   zeroQ    false if Q was computed;
 *   noqr     updated flag if ||dx_0|| is zero for the first time. */

void CPU_cmplx4_qrbs_head
 ( int dim, int degp1,
   double ***matrehihi, double ***matrelohi,
   double ***matrehilo, double ***matrelolo,
   double ***matimhihi, double ***matimlohi,
   double ***matimhilo, double ***matimlolo,
   double **rhsrehihi, double **rhsrelohi,
   double **rhsrehilo, double **rhsrelolo,
   double **rhsimhihi, double **rhsimlohi,
   double **rhsimhilo, double **rhsimlolo,
   double **solrehihi, double **solrelohi,
   double **solrehilo, double **solrelolo,
   double **solimhihi, double **solimlohi,
   double **solimhilo, double **solimlolo,
   double **Qrehihi, double **Qrelohi, double **Qrehilo, double **Qrelolo,
   double **Qimhihi, double **Qimlohi, double **Qimhilo, double **Qimlolo,
   double **Rrehihi, double **Rrelohi, double **Rrehilo, double **Rrelolo,
   double **Rimhihi, double **Rimlohi, double **Rimhilo, double **Rimlolo,
   double *wrkvecrehihi, double *wrkvecrelohi,
   double *wrkvecrehilo, double *wrkvecrelolo,
   double *wrkvecimhihi, double *wrkvecimlohi,
   double *wrkvecimhilo, double *wrkvecimlolo,
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
 *   matrehihi are degp1 matrices of dimension dim;
 *   matrelohi are degp1 matrices of dimension dim;
 *   matrehilo are degp1 matrices of dimension dim;
 *   matrelolo are degp1 matrices of dimension dim;
 *   matimhihi are degp1 matrices of dimension dim;
 *   matimlohi are degp1 matrices of dimension dim;
 *   matimhilo are degp1 matrices of dimension dim;
 *   matimlolo are degp1 matrices of dimension dim;
 *   rhsrehihi are degp1 vectors of dimension dim;
 *   rhsrelohi are degp1 vectors of dimension dim;
 *   rhsrehilo are degp1 vectors of dimension dim;
 *   rhsrelolo are degp1 vectors of dimension dim;
 *   rhsimhihi are degp1 vectors of dimension dim;
 *   rhsimlohi are degp1 vectors of dimension dim;
 *   rhsimhilo are degp1 vectors of dimension dim;
 *   rhsimlolo are degp1 vectors of dimension dim;
 *   solrehihi has space allocated for degp1 vectors of dimension dim;
 *   solrelohi has space allocated for degp1 vectors of dimension dim;
 *   solrehilo has space allocated for degp1 vectors of dimension dim;
 *   solrelolo has space allocated for degp1 vectors of dimension dim;
 *   solimhihi has space allocated for degp1 vectors of dimension dim;
 *   solimlohi has space allocated for degp1 vectors of dimension dim;
 *   solimhilo has space allocated for degp1 vectors of dimension dim;
 *   solimlolo has space allocated for degp1 vectors of dimension dim;
 *   Qrehihi  space allocated for a matrix of dimension dim;
 *   Qrelohi  space allocated for a matrix of dimension dim;
 *   Qrehilo  space allocated for a matrix of dimension dim;
 *   Qrelolo  space allocated for a matrix of dimension dim;
 *   Qimhihi  space allocated for a matrix of dimension dim;
 *   Qimlohi  space allocated for a matrix of dimension dim;
 *   Qimhilo  space allocated for a matrix of dimension dim;
 *   Qimlolo  space allocated for a matrix of dimension dim;
 *   Rrehihi  space allocated for a matrix of dimension dim;
 *   Rrelohi  space allocated for a matrix of dimension dim;
 *   Rrehilo  space allocated for a matrix of dimension dim;
 *   Rrelolo  space allocated for a matrix of dimension dim;
 *   Rimhihi  space allocated for a matrix of dimension dim;
 *   Rimlohi  space allocated for a matrix of dimension dim;
 *   Rimhilo  space allocated for a matrix of dimension dim;
 *   Rimlolo  space allocated for a matrix of dimension dim;
 *   wrkvecrehihi is work space allocated for a vector of dimension dim;
 *   wrkvecrelohi is work space allocated for a vector of dimension dim;
 *   wrkvecrehilo is work space allocated for a vector of dimension dim;
 *   wrkvecrelolo is work space allocated for a vector of dimension dim;
 *   wrkvecimhihi is work space allocated for a vector of dimension dim;
 *   wrkvecimlohi is work space allocated for a vector of dimension dim;
 *   wrkvecimhilo is work space allocated for a vector of dimension dim;
 *   wrkvecimlolo is work space allocated for a vector of dimension dim;
 *   zeroQ    if true, then Q is zero and Q must be computed;
 *   noqr     flag if true, then no qr, only when not zeroQ;
 *   vrblvl   the verbose level (0 for silent).
 *
 * ON RETURN :
 *   Qrehihi  highest doubles of the real parts of the Q in the QR;
 *   Qrelohi  second highest doubles of the real parts of the Q in the QR;
 *   Qrehilo  second lowest doubles of the real parts of the Q in the QR;
 *   Qrelolo  lowest doubles of the real parts of the Q in the QR;
 *   Qimhihi  highest doubles of the imaginary parts of the Q in the QR;
 *   Qimlohi  second highest doubles of the imaginary parts of the Q in the QR;
 *   Qimhilo  second lowest doubles of the imaginary parts of the Q in the QR;
 *   Qimlolo  lowest doubles of the imaginary parts of the Q in the QR;
 *   Rrehihi  highest doubles of the real parts of the R in the QR;
 *   Rrelohi  second highest doubles of the real parts of the R in the QR;
 *   Rrehilo  second lowest doubles of the real parts of the R in the QR;
 *   Rrelolo  lowest doubles of the real parts of the R in the QR;
 *   Rimhihi  highest doubles of the imaginary parts of the R in the QR;
 *   Rimlohi  second highest doubles of the imaginary parts of the R in the QR;
 *   Rimhilo  second lowest doubles of the imaginary parts of the R in the QR;
 *   Rimlolo  lowest doubles of the imaginary parts of the R in the QR;
 *   wrkvecrehihi is work space used to solve the linear system;
 *   wrkvecrelohi is work space used to solve the linear system;
 *   wrkvecrehilo is work space used to solve the linear system;
 *   wrkvecrelolo is work space used to solve the linear system;
 *   wrkvecimhihi is work space used to solve the linear system;
 *   wrkvecimlohi is work space used to solve the linear system;
 *   wrkvecimhilo is work space used to solve the linear system;
 *   wrkvecimlolo is work space used to solve the linear system;
 *   solrehihi are the highest doubles of the real parts of the head term;
 *   solrelohi are the second highest doubles of the real parts of the head;
 *   solrehilo are the second lowest doubles of the real parts of the head;
 *   solrelolo are the lowest doubles of the real parts of the head term;
 *   solimhihi are the highest doubles of the imag parts of the head term;
 *   solimlohi are the second highest doubles of the imag parts of the head;
 *   solimhilo are the second lowest doubles of the imag parts of the head;
 *   solimlolo are the lowest doubles of the imag parts of the head term;
 *   zeroQ    false if Q was computed;
 *   noqr     updated flag if ||dx_0|| is zero for the first time. */

void CPU_dbl4_qrbs_tail
 ( int dim, int degp1, int tailidx,
   double ***mathihi, double ***matlohi, double ***mathilo, double ***matlolo,
   double **rhshihi, double **rhslohi, double **rhshilo, double **rhslolo,
   double **solhihi, double **sollohi, double **solhilo, double **sollolo,
   double **Qhihi, double **Qlohi, double **Qhilo, double **Qlolo,
   double **Rhihi, double **Rlohi, double **Rhilo, double **Rlolo,
   double *wrkvechihi, double *wrkveclohi,
   double *wrkvechilo, double *wrkveclolo,
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
 *   mathihi  degp1 matrices of dimension dim;
 *   matlohi  degp1 matrices of dimension dim;
 *   mathilo  degp1 matrices of dimension dim;
 *   matlolo  degp1 matrices of dimension dim;
 *   rhshihi  degp1 vectors of dimension dim;
 *   rhslohi  degp1 vectors of dimension dim;
 *   rhshilo  degp1 vectors of dimension dim;
 *   rhslolo  degp1 vectors of dimension dim;
 *   solhihi  space allocated for degp1 vectors of dimension dim,
 *            with the leading highest doubles defined;
 *   sollohi  space allocated for degp1 vectors of dimension dim,
 *            with the leading second highest doubles defined;
 *   solhilo  space allocated for degp1 vectors of dimension dim,
 *            with the leading second lowest doubles defined;
 *   sollolo  space allocated for degp1 vectors of dimension dim,
 *            with the leading lowest doubles defined;
 *   Qhihi    highest doubles of the Q of the QR factorization;
 *   Qlohi    second highest doubles of the Q of the QR factorization;
 *   Qhilo    second lowest doubles of the Q of the QR factorization;
 *   Qlolo    lowest doubles of the Q of the QR factorization;
 *   Rhihi    highest doubles of the R of the QR factorization;
 *   Rlohi    second highest doubles of the R of the QR factorization;
 *   Rhilo    second lowest doubles of the R of the QR factorization;
 *   Rlolo    lowest doubles of the R of the QR factorization;
 *   wrkvechihi is work space vector of dimension dim for the substitution;
 *   wrkveclohi is work space vector of dimension dim for the substitution;
 *   wrkvechilo is work space vector of dimension dim for the substitution;
 *   wrkveclolo is work space vector of dimension dim for the substitution;
 *   vrblvl   the verbose level (0 for silent).
 *
 * ON RETURN :
 *   rhshihi  updated highest doubles of right hand side used as work space;
 *   rhslohi  updated second highest doubles of right hand side;
 *   rhshilo  updated second lowest doubles of right hand side;
 *   rhslolo  updated lowest doubles of right hand side used as work space;
 *   solhihi  highest doubles of all coefficients of the solution;
 *   sollohi  second highest doubles of all coefficients of the solution;
 *   solhilo  second lowest doubles of all coefficients of the solution;
 *   sollolo  lowest doubles of all coefficients of the solution;
 *   upidx    counts the number of updates skipped;
 *   bsidx    counts the number of backsubstitutions skipped;
 *   newtail  the new value for tailidx. */

void CPU_cmplx4_qrbs_tail
 ( int dim, int degp1, int tailidx,
   double ***matrehihi, double ***matrelohi,
   double ***matrehilo, double ***matrelolo,
   double ***matimhihi, double ***matimlohi,
   double ***matimhilo, double ***matimlolo,
   double **rhsrehihi, double **rhsrelohi,
   double **rhsrehilo, double **rhsrelolo,
   double **rhsimhihi, double **rhsimlohi,
   double **rhsimhilo, double **rhsimlolo,
   double **solrehihi, double **solrelohi,
   double **solrehilo, double **solrelolo,
   double **solimhihi, double **solimlohi,
   double **solimhilo, double **solimlolo,
   double **Qrehihi, double **Qrelohi, double **Qrehilo, double **Qrelolo,
   double **Qimhihi, double **Qimlohi, double **Qimhilo, double **Qimlolo,
   double **Rrehihi, double **Rrelohi, double **Rrehilo, double **Rrelolo, 
   double **Rimhihi, double **Rimlohi, double **Rimhilo, double **Rimlolo,
   double *wrkvecrehihi, double *wrkvecrelohi,
   double *wrkvecrehilo, double *wrkvecrelolo,
   double *wrkvecimhihi, double *wrkvecimlohi,
   double *wrkvecimhilo, double *wrkvecimlolo,
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
 *   matrehihi are degp1 matrices of dimension dim;
 *   matrelohi are degp1 matrices of dimension dim;
 *   matrehilo are degp1 matrices of dimension dim;
 *   matrelolo are degp1 matrices of dimension dim;
 *   matimhihi are degp1 matrices of dimension dim;
 *   matimlohi are degp1 matrices of dimension dim;
 *   matimhilo are degp1 matrices of dimension dim;
 *   matimlolo are degp1 matrices of dimension dim;
 *   rhsrehihi are degp1 vectors of dimension dim;
 *   rhsrelohi are degp1 vectors of dimension dim;
 *   rhsrehilo are degp1 vectors of dimension dim;
 *   rhsrelolo are degp1 vectors of dimension dim;
 *   rhsimhihi are degp1 vectors of dimension dim;
 *   rhsimlohi are degp1 vectors of dimension dim;
 *   rhsimhilo are degp1 vectors of dimension dim;
 *   rhsimlolo are degp1 vectors of dimension dim;
 *   solrehihi has space for degp1 vectors of dimension dim, with the
 *            highest doubles of the real parts of the head;
 *   solrelohi has space for degp1 vectors of dimension dim, with the
 *            second highest doubles of the real parts of the head;
 *   solrehilo has space for degp1 vectors of dimension dim, with the
 *            second lowest doubles of the real parts of the;
 *   solrelolo has space for degp1 vectors of dimension dim, with the
 *            lowest doubles of the real parts of the;
 *   solimhihi has space for degp1 vectors of dimension dim, with the
 *            highest doubles of the imaginary parts of the head;
 *   solimlohi has space for degp1 vectors of dimension dim, with the
 *            second highest doubles of the imaginary parts of the head;
 *   solimhilo has space for degp1 vectors of dimension dim, with the
 *            second lowest doubles of the imaginary parts of the head;
 *   solimlolo has space for degp1 vectors of dimension dim, with the
 *            lowest doubles of the imaginary parts of the head;
 *   Qrehihi  highest doubles of the real parts of the Q of the QR;
 *   Qrelohi  second highest doubles of the real parts of the Q of the QR;
 *   Qrehilo  second lowest doubles of the real parts of the Q of the QR;
 *   Qrelolo  lowest doubles of the real parts of the Q of the QR;
 *   Qimhihi  highest doubles of the imaginary parts of the Q of the QR;
 *   Qimlohi  second highest doubles of the imaginary parts of the Q of the QR;
 *   Qimhilo  second lowest doubles of the imaginary parts of the Q of the QR;
 *   Qimlolo  lowest doubles of the imaginary parts of the Q of the QR;
 *   Rrehihi  highest doubles of the real parts of the R of the QR;
 *   Rrelohi  second highest doubles of the real parts of the R of the QR;
 *   Rrehilo  second lowest doubles of the real parts of the R of the QR;
 *   Rrelolo  lowest doubles of the real parts of the R of the QR;
 *   Rimhihi  highest doubles of the imaginary parts of the R of the QR;
 *   Rimlohi  second highest doubles of the imaginary parts of the R of the QR;
 *   Rimhilo  second lowest doubles of the imaginary parts of the R of the QR;
 *   Rimlolo  lowest doubles of the imaginary parts of the R of the QR;
 *   wrkvecrehihi is work space of dimension dim for the substitution;
 *   wrkvecrelohi is work space of dimension dim for the substitution;
 *   wrkvecrehilo is work space of dimension dim for the substitution;
 *   wrkvecrelolo is work space of dimension dim for the substitution;
 *   wrkvecimhihi is work space of dimension dim for the substitution;
 *   wrkvecimlohi is work space of dimension dim for the substitution;
 *   wrkvecimhilo is work space of dimension dim for the substitution;
 *   wrkvecimlolo is work space of dimension dim for the substitution;
 *   vrblvl   the verbose level (0 for silent).
 *
 * ON RETURN :
 *   rhsrehihi are the highest doubles of the real parts
 *            of the updated right hand side;
 *   rhsrelohi are the second highest doubles of the real parts
 *            of the updated right hand side;
 *   rhsimhilo are the second lowest doubles of the imaginary parts
 *            of the updated right hand side;
 *   rhsimlolo are the lowest doubles of the imaginary parts
 *            of the updated right hand side;
 *   solrehihi are the highest doubles of the real parts of the solution;
 *   solrelohi are the 2nd highest doubles of the real parts of the solution;
 *   solrehilo are the 2nd lowest doubles of the real parts of the solution;
 *   solrelolo are the lowest doubles of the real parts of the solution;
 *   solimhihi are the highest doubles of the imaginary parts of the solution;
 *   solimlohi are the 2nd highest doubles of the imag parts of the solution;
 *   solimhilo are the 2nd lowest doubles of the imag parts of the solution;
 *   solimlolo are the lowest doubles of the imaginary parts of the solution;
 *   upidx    counts the number of updates skipped;
 *   bsidx    counts the number of backsubstitutions skipped;
 *   newtail  the new value for tailidx. */

void CPU_dbl4_qrbs_solve
 ( int dim, int degp1, int tailidx,
   double ***mathihi, double ***matlohi, double ***mathilo, double ***matlolo,
   double **rhshihi, double **rhslohi, double **rhshilo, double **rhslolo,
   double **solhihi, double **sollohi, double **solhilo, double **sollolo,
   double **Qhihi, double **Qlohi, double **Qhilo, double **Qlolo,
   double **Rhihi, double **Rlohi, double **Rhilo, double **Rlolo,
   double *wrkvechihi, double *wrkveclohi,
   double *wrkvechilo, double *wrkveclolo,
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
 *   mathihi  degp1 matrices of dimension dim;
 *   matlohi  degp1 matrices of dimension dim;
 *   mathilo  degp1 matrices of dimension dim;
 *   matlolo  degp1 matrices of dimension dim;
 *   rhshihi  degp1 vectors of dimension dim;
 *   rhslohi  degp1 vectors of dimension dim;
 *   rhshilo  degp1 vectors of dimension dim;
 *   rhslolo  degp1 vectors of dimension dim;
 *   solhihi  space allocated for degp1 vectors of dimension dim;
 *   sollohi  space allocated for degp1 vectors of dimension dim;
 *   solhilo  space allocated for degp1 vectors of dimension dim;
 *   sollolo  space allocated for degp1 vectors of dimension dim;
 *   Qhihi    space allocated for a matrix of dimension dim;
 *   Qlohi    space allocated for a matrix of dimension dim;
 *   Qhilo    space allocated for a matrix of dimension dim;
 *   Qlolo    space allocated for a matrix of dimension dim;
 *   Rhihi    space allocated for a matrix of dimension dim;
 *   Rlohi    space allocated for a matrix of dimension dim;
 *   Rhilo    space allocated for a matrix of dimension dim;
 *   Rlolo    space allocated for a matrix of dimension dim;
 *   wrkvechihi has work space allocated for a vector of dimension dim;
 *   wrkveclohi has work space allocated for a vector of dimension dim;
 *   wrkvechilo has work space allocated for a vector of dimension dim;
 *   wrkveclolo has work space allocated for a vector of dimension dim;
 *   zeroQ    if true, then Q is zero and Q must be computed;
 *   noqr     flag if true, then no qr, only when not zeroQ;
 *   vrblvl   the verbose level (0 for silent).
 *
 * ON RETURN :
 *   Qhihi    highest doubles of the Q of the QR of the Jacobian matrix;
 *   Qlohi    second highest doubles of the Q of the QR of the Jacobian matrix;
 *   Qhilo    second lowest doubles of the Q of the QR of the Jacobian matrix;
 *   Qlolo    lowest doubles of the Q of the QR of the Jacobian matrix;
 *   Rhihi    highest doubles of the R in the QR of the Jacobian matrix;
 *   Rlohi    second highest doubles of the R in the QR of the Jacobian matrix;
 *   Rhilo    second lowest doubles of the R in the QR of the Jacobian matrix;
 *   Rlolo    lowest doubles of the R in the QR of the Jacobian matrix;
 *   wrkvechihi is work space used to solve the linear system;
 *   wrkveclohi is work space used to solve the linear system;
 *   wrkvechilo is work space used to solve the linear system;
 *   wrkveclolo is work space used to solve the linear system;
 *   solhihi  highest double coefficients of the solution series;
 *   sollohi  second lowest double coefficients of the solution series;
 *   solhilo  second highest double coefficients of the solution series;
 *   sollolo  lowest double coefficients of the solution series;
 *   zeroQ    false if Q was computed;
 *   noqr     updated flag if ||dx_0|| is zero for the first time;
 *   upidx    counts the number of updates skipped;
 *   bsidx    counts the number of backsubstitutions skipped;
 *   newtail  the new value for tailidx. */

void CPU_cmplx4_qrbs_solve
 ( int dim, int degp1, int tailidx,
   double ***matrehihi, double ***matrelohi,
   double ***matrehilo, double ***matrelolo,
   double ***matimhihi, double ***matimlohi, 
   double ***matimhilo, double ***matimlolo, 
   double **rhsrehihi, double **rhsrelohi,
   double **rhsrehilo, double **rhsrelolo,
   double **rhsimhihi, double **rhsimlohi,
   double **rhsimhilo, double **rhsimlolo,
   double **solrehihi, double **solrelohi,
   double **solrehilo, double **solrelolo,
   double **solimhihi, double **solimlohi,
   double **solimhilo, double **solimlolo,
   double **Qrehihi, double **Qrelohi, double **Qrehilo, double **Qrelolo,
   double **Qimhihi, double **Qimlohi, double **Qimhilo, double **Qimlolo,
   double **Rrehihi, double **Rrelohi, double **Rrehilo, double **Rrelolo,
   double **Rimhihi, double **Rimlohi, double **Rimhilo, double **Rimlolo,
   double *wrkvecrehihi, double *wrkvecrelohi,
   double *wrkvecrehilo, double *wrkvecrelolo,
   double *wrkvecimhihi, double *wrkvecimlohi,
   double *wrkvecimhilo, double *wrkvecimlolo,
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
 *   matrehihi are degp1 matrices of dimension dim;
 *   matrelohi are degp1 matrices of dimension dim;
 *   matrehilo are degp1 matrices of dimension dim;
 *   matrelolo are degp1 matrices of dimension dim;
 *   matimhihi are degp1 matrices of dimension dim;
 *   matimlohi are degp1 matrices of dimension dim;
 *   matimhilo are degp1 matrices of dimension dim;
 *   matimlolo are degp1 matrices of dimension dim;
 *   rhsrehihi are degp1 vectors of dimension dim;
 *   rhsrelohi are degp1 vectors of dimension dim;
 *   rhsrehilo are degp1 vectors of dimension dim;
 *   rhsrelolo are degp1 vectors of dimension dim;
 *   rhsimhihi are degp1 vectors of dimension dim;
 *   rhsimlohi are degp1 vectors of dimension dim;
 *   rhsimhilo are degp1 vectors of dimension dim;
 *   rhsimlolo are degp1 vectors of dimension dim;
 *   solrehihi has space allocated for degp1 vectors of dimension dim;
 *   solrelohi has space allocated for degp1 vectors of dimension dim;
 *   solrehilo has space allocated for degp1 vectors of dimension dim;
 *   solrelolo has space allocated for degp1 vectors of dimension dim;
 *   solimhihi has space allocated for degp1 vectors of dimension dim;
 *   solimlohi has space allocated for degp1 vectors of dimension dim;
 *   solimhilo has space allocated for degp1 vectors of dimension dim;
 *   solimlolo has space allocated for degp1 vectors of dimension dim;
 *   Qrehihi  space allocated for a matrix of dimension dim;
 *   Qrelohi  space allocated for a matrix of dimension dim;
 *   Qrehilo  space allocated for a matrix of dimension dim;
 *   Qrelolo  space allocated for a matrix of dimension dim;
 *   Qimhihi  space allocated for a matrix of dimension dim;
 *   Qimlohi  space allocated for a matrix of dimension dim;
 *   Qimhilo  space allocated for a matrix of dimension dim;
 *   Qimlolo  space allocated for a matrix of dimension dim;
 *   Rrehihi  space allocated for a matrix of dimension dim;
 *   Rrelohi  space allocated for a matrix of dimension dim;
 *   Rrehilo  space allocated for a matrix of dimension dim;
 *   Rrelolo  space allocated for a matrix of dimension dim;
 *   Rimhihi  space allocated for a matrix of dimension dim;
 *   Rimlohi  space allocated for a matrix of dimension dim;
 *   Rimhilo  space allocated for a matrix of dimension dim;
 *   Rimlolo  space allocated for a matrix of dimension dim;
 *   wrkvecrehihi has work space allocated for a vector of dimension dim;
 *   wrkvecrelohi has work space allocated for a vector of dimension dim;
 *   wrkvecrehilo has work space allocated for a vector of dimension dim;
 *   wrkvecrelolo has work space allocated for a vector of dimension dim;
 *   wrkvecimhihi has work space allocated for a vector of dimension dim;
 *   wrkvecimlohi has work space allocated for a vector of dimension dim;
 *   wrkvecimhilo has work space allocated for a vector of dimension dim;
 *   wrkvecimlolo has work space allocated for a vector of dimension dim;
 *   zeroQ    if true, then Q is zero and Q must be computed;
 *   noqr     flag if true, then no qr, only when not zeroQ;
 *   vrblvl   the verbose level (0 for silent).
 *
 * ON RETURN :
 *   Qrehihi  highest doubles of the real parts of the Q in the QR;
 *   Qrelohi  second highest doubles of the real parts of the Q in the QR;
 *   Qrehilo  second lowest doubles of the real parts of the Q in the QR;
 *   Qrelolo  lowest doubles of the real parts of the Q in the QR;
 *   Qimhihi  highest doubles of the imaginary parts of the Q in the QR;
 *   Qimlohi  second highest doubles of the imaginary parts of the Q in the QR;
 *   Qimhilo  second lowest doubles of the imaginary parts of the Q in the QR;
 *   Qimlolo  lowest doubles of the imaginary parts of the Q in the QR;
 *   Rrehihi  highest doubles of the real parts of the R in the QR;
 *   Rrelohi  second highest doubles of the real parts of the R in the QR;
 *   Rrehilo  second lowest doubles of the real parts of the R in the QR;
 *   Rrelolo  lowest doubles of the real parts of the R in the QR;
 *   Rimhihi  highest doubles of the imaginary parts of the R in the QR;
 *   Rimlohi  second highest doubles of the imaginary parts of the R in the QR;
 *   Rimhilo  second lowest doubles of the imaginary parts of the R in the QR;
 *   Rimlolo  lowest doubles of the imaginary parts of the R in the QR;
 *   wrkvecrehihi is work space used to solve the linear systems;
 *   wrkvecrelohi is work space used to solve the linear systems;
 *   wrkvecrehilo is work space used to solve the linear systems;
 *   wrkvecrelolo is work space used to solve the linear systems;
 *   wrkvecimhihi is work space used to solve the linear systems;
 *   wrkvecimlohi is work space used to solve the linear systems;
 *   wrkvecimhilo is work space used to solve the linear systems;
 *   wrkvecimlolo is work space used to solve the linear systems;
 *   solrehihi are the highest doubles of the real parts of the solution;
 *   solrelohi are the second highest doubles of the real parts of the solution;
 *   solrehilo are the second lowest doubles of the real parts of the solution;
 *   solrelolo are the lowest doubles of the real parts of the solution;
 *   solimhihi are the highest doubles of the imag parts of the solution;
 *   solimlohi are the second highest doubles of the imag parts of the solution;
 *   solimhilo are the second lowest doubles of the imag parts of the solution;
 *   solimlolo are the lowest doubles of the imag parts of the solution;
 *   zeroQ    false if Q was computed;
 *   noqr     updated flag if ||dx_0|| is zero for the first time;
 *   upidx    counts the number of updates skipped;
 *   bsidx    counts the number of backsubstitutions skipped;
 *   newtail  the new value for tailidx. */

void CPU_dbl4_linear_residue
 ( int dim, int degp1, int tailidx,
   double ***mathihi, double ***matlohi, double ***mathilo, double ***matlolo, 
   double **rhshihi, double **rhslohi, double **rhshilo, double **rhslolo,
   double **solhihi, double **sollohi, double **solhilo, double **sollolo,
   double **resvechihi, double **resveclohi,
   double **resvechilo, double **resveclolo,
   double *resmaxhihi, double *resmaxlohi,
   double *resmaxhilo, double *resmaxlolo, int vrblvl );
/*
 * DESCRIPTION :
 *   Computes the residual of the linear power series system.
 *
 * ON ENTRY :
 *   dim      the dimension of the matrices and vectors;
 *   degp1    degree plus one, the size of the matrix system;
 *   tailidx  the index of the start of the update in the tail;
 *   mathihi  degp1 matrices of dimension dim;
 *   matlohi  degp1 matrices of dimension dim;
 *   mathilo  degp1 matrices of dimension dim;
 *   matlolo  degp1 matrices of dimension dim;
 *   rhshihi  degp1 right hand side vectors of dimension dim;
 *   rhslohi  degp1 right hand side vectors of dimension dim;
 *   rhshilo  degp1 right hand side vectors of dimension dim;
 *   rhslolo  degp1 right hand side vectors of dimension dim;
 *   solhihi  degp1 solution vectors of dimension dim;
 *   sollohi  degp1 solution vectors of dimension dim;
 *   solhilo  degp1 solution vectors of dimension dim;
 *   sollolo  degp1 solution vectors of dimension dim;
 *   resvechihi has space for the residual power series;
 *   resveclohi has space for the residual power series;
 *   resvechilo has space for the residual power series;
 *   resveclolo has space for the residual power series;
 *   vrblvl   is the verbose level.
 *
 * ON RETURN :
 *   resvechihi are the highest doubles of the residual power series;
 *   resveclohi are the 2nd highest doubles of the residual power series;
 *   resvechilo are the 2nd lowest doubles of the residual power series;
 *   resveclolo are the lowest doubles of the residual power series;
 *   resmaxhihi is the highest double of the maximum component in resvec;
 *   resmaxlohi is the 2nd highest double of the maximum component in resvec;
 *   resmaxhilo is the 2nd lowest double of the maximum component in resvec;
 *   resmaxlolo is the lowest double of the maximum component in resvec. */

void CPU_cmplx4_linear_residue
 ( int dim, int degp1, int tailidx,
   double ***matrehihi, double ***matrelohi,
   double ***matrehilo, double ***matrelolo,
   double ***matimhihi, double ***matimlohi,
   double ***matimhilo, double ***matimlolo,
   double **rhsrehihi, double **rhsrelohi,
   double **rhsrehilo, double **rhsrelolo,
   double **rhsimhihi, double **rhsimlohi, 
   double **rhsimhilo, double **rhsimlolo, 
   double **solrehihi, double **solrelohi,
   double **solrehilo, double **solrelolo,
   double **solimhihi, double **solimlohi,
   double **solimhilo, double **solimlolo,
   double **resvecrehihi, double **resvecrelohi,
   double **resvecrehilo, double **resvecrelolo,
   double **resvecimhihi, double **resvecimlohi,
   double **resvecimhilo, double **resvecimlolo,
   double *resmaxhihi, double *resmaxlohi,
   double *resmaxhilo, double *resmaxlolo, int vrblvl );
/*
 * DESCRIPTION :
 *   Computes the residual of the linear power series system.
 *
 * ON ENTRY :
 *   dim      the dimension of the matrices and vectors;
 *   degp1    degree plus one, the size of the matrix system;
 *   tailidx  the index of the start of the update in the tail;
 *   matrehihi are degp1 matrices of dimension dim;
 *   matrelohi are degp1 matrices of dimension dim;
 *   matrehilo are degp1 matrices of dimension dim;
 *   matrelolo are degp1 matrices of dimension dim;
 *   matimhihi are degp1 matrices of dimension dim;
 *   matimlohi are degp1 matrices of dimension dim;
 *   matimhilo are degp1 matrices of dimension dim;
 *   matimlolo are degp1 matrices of dimension dim;
 *   rhsrehihi are degp1 right hand side vectors of dimension dim;
 *   rhsrelohi are degp1 right hand side vectors of dimension dim;
 *   rhsrehilo are degp1 right hand side vectors of dimension dim;
 *   rhsrelolo are degp1 right hand side vectors of dimension dim;
 *   rhsimhihi are degp1 right hand side vectors of dimension dim;
 *   rhsimlohi are degp1 right hand side vectors of dimension dim;
 *   rhsimhilo are degp1 right hand side vectors of dimension dim;
 *   rhsimlolo are degp1 right hand side vectors of dimension dim;
 *   solrehihi are degp1 solution vectors of dimension dim;
 *   solrelohi are degp1 solution vectors of dimension dim;
 *   solrehilo are degp1 solution vectors of dimension dim;
 *   solrelolo are degp1 solution vectors of dimension dim;
 *   solimhihi are degp1 solution vectors of dimension dim;
 *   solimlohi are degp1 solution vectors of dimension dim;
 *   solimhilo are degp1 solution vectors of dimension dim;
 *   solimlolo are degp1 solution vectors of dimension dim;
 *   resvecrehihi has space for the residual power series;
 *   resvecrelohi has space for the residual power series;
 *   resvecrehilo has space for the residual power series;
 *   resvecrelolo has space for the residual power series;
 *   resvecimhihi has space for the residual power series;
 *   resvecimlohi has space for the residual power series;
 *   resvecimhilo has space for the residual power series;
 *   resvecimlolo has space for the residual power series;
 *   vrblvl   is the verbose level.
 *
 * ON RETURN :
 *   resvecrehihi are the highest doubles of the real parts
 *            of the residual series;
 *   resvecrelohi are the second highest doubles of the real parts
 *            of the residual series;
 *   resvecrehilo are the second lowest doubles of the real parts
 *            of the residual series;
 *   resvecrelolo are the lowest doubles of the real parts
 *            of the residual series;
 *   resvecimhihi are the highest doubles of the imaginary parts
 *            of the residual series;
 *   resvecimlohi are the second highest doubles of the imaginary parts
 *            of the residual series;
 *   resvecimhilo are the second lowest doubles of the imaginary parts
 *            of the residual series;
 *   resvecimlolo are the lowest doubles of the imaginary parts
 *            of the residual series;
 *   resmaxhihi is the highest double of the max norm of the residual;
 *   resmaxlohi is the second highest double of the max norm of the residual;
 *   resmaxhilo is the second lowest double of the max norm of the residual;
 *   resmaxlolo is the lowest double of the max norm of the residual. */

#endif
