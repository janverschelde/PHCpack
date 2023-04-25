// The file dbl8_bals_host.h specifies functions to solve linear systems
// of power series by linearization in octo double precision.

#ifndef __dbl8_bals_host_h__
#define __dbl8_bals_host_h__

void CPU_dbl8_qrbs_head
 ( int dim, int degp1,
   double ***mathihihi, double ***matlohihi,
   double ***mathilohi, double ***matlolohi,
   double ***mathihilo, double ***matlohilo,
   double ***mathilolo, double ***matlololo,
   double **rhshihihi, double **rhslohihi,
   double **rhshilohi, double **rhslolohi,
   double **rhshihilo, double **rhslohilo,
   double **rhshilolo, double **rhslololo,
   double **solhihihi, double **sollohihi,
   double **solhilohi, double **sollolohi,
   double **solhihilo, double **sollohilo,
   double **solhilolo, double **sollololo,
   double **Qhihihi, double **Qlohihi, double **Qhilohi, double **Qlolohi,
   double **Qhihilo, double **Qlohilo, double **Qhilolo, double **Qlololo,
   double **Rhihihi, double **Rlohihi, double **Rhilohi, double **Rlolohi, 
   double **Rhihilo, double **Rlohilo, double **Rhilolo, double **Rlololo, 
   double *wrkvechihihi, double *wrkveclohihi,
   double *wrkvechilohi, double *wrkveclolohi,
   double *wrkvechihilo, double *wrkveclohilo,
   double *wrkvechilolo, double *wrkveclololo,
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
 *   mathihihi are degp1 matrices of dimension dim;
 *   matlohihi are degp1 matrices of dimension dim;
 *   mathilohi are degp1 matrices of dimension dim;
 *   matlolohi are degp1 matrices of dimension dim;
 *   mathihilo are degp1 matrices of dimension dim;
 *   matlohilo are degp1 matrices of dimension dim;
 *   mathilolo are degp1 matrices of dimension dim;
 *   matlololo are degp1 matrices of dimension dim;
 *   rhshihihi are degp1 vectors of dimension dim;
 *   rhslohihi are degp1 vectors of dimension dim;
 *   rhshilohi are degp1 vectors of dimension dim;
 *   rhslolohi are degp1 vectors of dimension dim;
 *   rhshihilo are degp1 vectors of dimension dim;
 *   rhslohilo are degp1 vectors of dimension dim;
 *   rhshilolo are degp1 vectors of dimension dim;
 *   rhslololo are degp1 vectors of dimension dim;
 *   solhihihi has space allocated for degp1 vectors of dimension dim;
 *   sollohihi has space allocated for degp1 vectors of dimension dim;
 *   solhilohi has space allocated for degp1 vectors of dimension dim;
 *   sollolohi has space allocated for degp1 vectors of dimension dim;
 *   solhihilo has space allocated for degp1 vectors of dimension dim;
 *   sollohilo has space allocated for degp1 vectors of dimension dim;
 *   solhilolo has space allocated for degp1 vectors of dimension dim;
 *   sollololo has space allocated for degp1 vectors of dimension dim;
 *   wrkvechihihi is work space allocated for a vector of dimension dim;
 *   wrkveclohihi is work space allocated for a vector of dimension dim;
 *   wrkvechilohi is work space allocated for a vector of dimension dim;
 *   wrkveclolohi is work space allocated for a vector of dimension dim;
 *   wrkvechihilo is work space allocated for a vector of dimension dim;
 *   wrkveclohilo is work space allocated for a vector of dimension dim;
 *   wrkvechilolo is work space allocated for a vector of dimension dim;
 *   wrkveclololo is work space allocated for a vector of dimension dim;
 *   Qhihihi  space allocated for a matrix of dimension dim;
 *   Qlohihi  space allocated for a matrix of dimension dim;
 *   Qhilohi  space allocated for a matrix of dimension dim;
 *   Qlolohi  space allocated for a matrix of dimension dim;
 *   Qhihilo  space allocated for a matrix of dimension dim;
 *   Qlohilo  space allocated for a matrix of dimension dim;
 *   Qhilolo  space allocated for a matrix of dimension dim;
 *   Qlololo  space allocated for a matrix of dimension dim;
 *   Rhihihi  space allocated for a matrix of dimension dim;
 *   Rlohihi  space allocated for a matrix of dimension dim;
 *   Rhilohi  space allocated for a matrix of dimension dim;
 *   Rlolohi  space allocated for a matrix of dimension dim;
 *   Rhihilo  space allocated for a matrix of dimension dim;
 *   Rlohilo  space allocated for a matrix of dimension dim;
 *   Rhilolo  space allocated for a matrix of dimension dim;
 *   Rlololo  space allocated for a matrix of dimension dim;
 *   wrkvechihihi is work space allocated for a vector of dimension dim;
 *   wrkveclohihi is work space allocated for a vector of dimension dim;
 *   wrkvechilohi is work space allocated for a vector of dimension dim;
 *   wrkveclolohi is work space allocated for a vector of dimension dim;
 *   wrkvechihilo is work space allocated for a vector of dimension dim;
 *   wrkveclohilo is work space allocated for a vector of dimension dim;
 *   wrkvechilolo is work space allocated for a vector of dimension dim;
 *   wrkveclololo is work space allocated for a vector of dimension dim;
 *   zeroQ    if true, then Q is zero and Q must be computed;
 *   noqr     flag if true, then no qr, only when not zeroQ;
 *   vrblvl   the verbose level (0 for silent).
 *
 * ON RETURN :
 *   Qhihihi  highest doubles of the Q of the QR of the Jacobian;
 *   Qlohihi  second highest doubles of the Q of the QR of the Jacobian;
 *   Qhilohi  third highest doubles of the Q of the QR of the Jacobian;
 *   Qlolohi  fourth highest doubles of the Q of the QR of the Jacobian;
 *   Qhihilo  fourth lowest doubles of the Q of the QR of the Jacobian;
 *   Qlohilo  third lowest doubles of the Q of the QR of the Jacobian;
 *   Qhilolo  second lowest doubles of the Q of the QR of the Jacobian;
 *   Qlololo  lowest doubles of the Q of the QR of the Jacobian;
 *   Rhihihi  highest doubles of the R in the QR of the Jacobian;
 *   Rlohihi  second highest doubles of the R in the QR of the Jacobian;
 *   Rhilohi  third highest doubles of the R in the QR of the Jacobian;
 *   Rlolohi  fourth highest doubles of the R in the QR of the Jacobian;
 *   Rhihilo  fourth lowest doubles of the R in the QR of the Jacobian;
 *   Rlohilo  third lowest doubles of the R in the QR of the Jacobian;
 *   Rhilolo  second lowest doubles of the R in the QR of the Jacobian;
 *   Rlololo  lowest doubles of the R in the QR of the Jacobian;
 *   wrkvechihihi is work space used to solve the linear system;
 *   wrkveclohihi is work space used to solve the linear system;
 *   wrkvechilohi is work space used to solve the linear system;
 *   wrkveclolohi is work space used to solve the linear system;
 *   wrkvechihilo is work space used to solve the linear system;
 *   wrkveclohilo is work space used to solve the linear system;
 *   wrkvechilolo is work space used to solve the linear system;
 *   wrkveclololo is work space used to solve the linear system;
 *   solhihihi  highest doubles of the coefficients of the solution;
 *   sollohihi  second highest doubles of the coefficients of the solution;
 *   solhilohi  third highest doubles of the coefficients of the solution;
 *   sollolohi  fourth highest doubles of the coefficients of the solution;
 *   solhihilo  fourth lowest doubles of the coefficients of the solution;
 *   sollohilo  third lowest doubles of the coefficients of the solution;
 *   solhilolo  second lowest doubles of the coefficients of the solution;
 *   sollololo  lowest doubles of the coefficients of the solution;
 *   zeroQ    false if Q was computed;
 *   noqr     updated flag if ||dx_0|| is zero for the first time. */

void CPU_cmplx8_qrbs_head
 ( int dim, int degp1,
   double ***matrehihihi, double ***matrelohihi,
   double ***matrehilohi, double ***matrelolohi,
   double ***matrehihilo, double ***matrelohilo,
   double ***matrehilolo, double ***matrelololo,
   double ***matimhihihi, double ***matimlohihi,
   double ***matimhilohi, double ***matimlolohi,
   double ***matimhihilo, double ***matimlohilo,
   double ***matimhilolo, double ***matimlololo,
   double **rhsrehihihi, double **rhsrelohihi,
   double **rhsrehilohi, double **rhsrelolohi,
   double **rhsrehihilo, double **rhsrelohilo,
   double **rhsrehilolo, double **rhsrelololo,
   double **rhsimhihihi, double **rhsimlohihi,
   double **rhsimhilohi, double **rhsimlolohi,
   double **rhsimhihilo, double **rhsimlohilo,
   double **rhsimhilolo, double **rhsimlololo,
   double **solrehihihi, double **solrelohihi,
   double **solrehilohi, double **solrelolohi,
   double **solrehihilo, double **solrelohilo,
   double **solrehilolo, double **solrelololo,
   double **solimhihihi, double **solimlohihi,
   double **solimhilohi, double **solimlolohi,
   double **solimhihilo, double **solimlohilo,
   double **solimhilolo, double **solimlololo,
   double **Qrehihihi, double **Qrelohihi,
   double **Qrehilohi, double **Qrelolohi,
   double **Qrehihilo, double **Qrelohilo,
   double **Qrehilolo, double **Qrelololo,
   double **Qimhihihi, double **Qimlohihi,
   double **Qimhilohi, double **Qimlolohi,
   double **Qimhihilo, double **Qimlohilo,
   double **Qimhilolo, double **Qimlololo,
   double **Rrehihihi, double **Rrelohihi,
   double **Rrehilohi, double **Rrelolohi,
   double **Rrehihilo, double **Rrelohilo,
   double **Rrehilolo, double **Rrelololo,
   double **Rimhihihi, double **Rimlohihi,
   double **Rimhilohi, double **Rimlolohi,
   double **Rimhihilo, double **Rimlohilo,
   double **Rimhilolo, double **Rimlololo,
   double *wrkvecrehihihi, double *wrkvecrelohihi,
   double *wrkvecrehilohi, double *wrkvecrelolohi,
   double *wrkvecrehihilo, double *wrkvecrelohilo,
   double *wrkvecrehilolo, double *wrkvecrelololo,
   double *wrkvecimhihihi, double *wrkvecimlohihi,
   double *wrkvecimhilohi, double *wrkvecimlolohi,
   double *wrkvecimhihilo, double *wrkvecimlohilo,
   double *wrkvecimhilolo, double *wrkvecimlololo,
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
 *   matrehihihi are degp1 matrices of dimension dim;
 *   matrelohihi are degp1 matrices of dimension dim;
 *   matrehilohi are degp1 matrices of dimension dim;
 *   matrelolohi are degp1 matrices of dimension dim;
 *   matrehihilo are degp1 matrices of dimension dim;
 *   matrelohilo are degp1 matrices of dimension dim;
 *   matrehilolo are degp1 matrices of dimension dim;
 *   matrelololo are degp1 matrices of dimension dim;
 *   matimhihihi are degp1 matrices of dimension dim;
 *   matimlohihi are degp1 matrices of dimension dim;
 *   matimhilohi are degp1 matrices of dimension dim;
 *   matimlolohi are degp1 matrices of dimension dim;
 *   matimhihilo are degp1 matrices of dimension dim;
 *   matimlohilo are degp1 matrices of dimension dim;
 *   matimhilolo are degp1 matrices of dimension dim;
 *   matimlololo are degp1 matrices of dimension dim;
 *   rhsrehihihi are degp1 vectors of dimension dim;
 *   rhsrelohihi are degp1 vectors of dimension dim;
 *   rhsrehilohi are degp1 vectors of dimension dim;
 *   rhsrelolohi are degp1 vectors of dimension dim;
 *   rhsrehihilo are degp1 vectors of dimension dim;
 *   rhsrelohilo are degp1 vectors of dimension dim;
 *   rhsrehilolo are degp1 vectors of dimension dim;
 *   rhsrelololo are degp1 vectors of dimension dim;
 *   rhsimhihihi are degp1 vectors of dimension dim;
 *   rhsimlohihi are degp1 vectors of dimension dim;
 *   rhsimhilohi are degp1 vectors of dimension dim;
 *   rhsimlolohi are degp1 vectors of dimension dim;
 *   rhsimhihilo are degp1 vectors of dimension dim;
 *   rhsimlohilo are degp1 vectors of dimension dim;
 *   rhsimhilolo are degp1 vectors of dimension dim;
 *   rhsimlololo are degp1 vectors of dimension dim;
 *   solrehihihi has space allocated for degp1 vectors of dimension dim;
 *   solrelohihi has space allocated for degp1 vectors of dimension dim;
 *   solrehilohi has space allocated for degp1 vectors of dimension dim;
 *   solrelolohi has space allocated for degp1 vectors of dimension dim;
 *   solrehihilo has space allocated for degp1 vectors of dimension dim;
 *   solrelohilo has space allocated for degp1 vectors of dimension dim;
 *   solrehilolo has space allocated for degp1 vectors of dimension dim;
 *   solrelololo has space allocated for degp1 vectors of dimension dim;
 *   solimhihihi has space allocated for degp1 vectors of dimension dim;
 *   solimlohihi has space allocated for degp1 vectors of dimension dim;
 *   solimhilohi has space allocated for degp1 vectors of dimension dim;
 *   solimlolohi has space allocated for degp1 vectors of dimension dim;
 *   solimhihilo has space allocated for degp1 vectors of dimension dim;
 *   solimlohilo has space allocated for degp1 vectors of dimension dim;
 *   solimhilolo has space allocated for degp1 vectors of dimension dim;
 *   solimlololo has space allocated for degp1 vectors of dimension dim;
 *   Qrehihihi has space allocated for a matrix of dimension dim;
 *   Qrelohihi has space allocated for a matrix of dimension dim;
 *   Qrehilohi has space allocated for a matrix of dimension dim;
 *   Qrelolohi has space allocated for a matrix of dimension dim;
 *   Qrehihilo has space allocated for a matrix of dimension dim;
 *   Qrelohilo has space allocated for a matrix of dimension dim;
 *   Qrehilolo has space allocated for a matrix of dimension dim;
 *   Qrelololo has space allocated for a matrix of dimension dim;
 *   Qimhihihi has space allocated for a matrix of dimension dim;
 *   Qimlohihi has space allocated for a matrix of dimension dim;
 *   Qimhilohi has space allocated for a matrix of dimension dim;
 *   Qimlolohi has space allocated for a matrix of dimension dim;
 *   Qimhihilo has space allocated for a matrix of dimension dim;
 *   Qimlohilo has space allocated for a matrix of dimension dim;
 *   Qimhilolo has space allocated for a matrix of dimension dim;
 *   Qimlololo has space allocated for a matrix of dimension dim;
 *   Rrehihihi has space allocated for a matrix of dimension dim;
 *   Rrelohihi has space allocated for a matrix of dimension dim;
 *   Rrehilohi has space allocated for a matrix of dimension dim;
 *   Rrelolohi has space allocated for a matrix of dimension dim;
 *   Rrehihilo has space allocated for a matrix of dimension dim;
 *   Rrelohilo has space allocated for a matrix of dimension dim;
 *   Rrehilolo has space allocated for a matrix of dimension dim;
 *   Rrelololo has space allocated for a matrix of dimension dim;
 *   Rimhihihi has space allocated for a matrix of dimension dim;
 *   Rimlohihi has space allocated for a matrix of dimension dim;
 *   Rimhilohi has space allocated for a matrix of dimension dim;
 *   Rimlolohi has space allocated for a matrix of dimension dim;
 *   Rimhihilo has space allocated for a matrix of dimension dim;
 *   Rimlohilo has space allocated for a matrix of dimension dim;
 *   Rimhilolo has space allocated for a matrix of dimension dim;
 *   Rimlololo has space allocated for a matrix of dimension dim;
 *   wrkvecrehihihi is work space allocated for a vector of dimension dim;
 *   wrkvecrelohihi is work space allocated for a vector of dimension dim;
 *   wrkvecrehilohi is work space allocated for a vector of dimension dim;
 *   wrkvecrelolohi is work space allocated for a vector of dimension dim;
 *   wrkvecrehihilo is work space allocated for a vector of dimension dim;
 *   wrkvecrelohilo is work space allocated for a vector of dimension dim;
 *   wrkvecrehilolo is work space allocated for a vector of dimension dim;
 *   wrkvecrelololo is work space allocated for a vector of dimension dim;
 *   wrkvecimhihihi is work space allocated for a vector of dimension dim;
 *   wrkvecimlohihi is work space allocated for a vector of dimension dim;
 *   wrkvecimhilohi is work space allocated for a vector of dimension dim;
 *   wrkvecimlolohi is work space allocated for a vector of dimension dim;
 *   wrkvecimhihilo is work space allocated for a vector of dimension dim;
 *   wrkvecimlohilo is work space allocated for a vector of dimension dim;
 *   wrkvecimhilolo is work space allocated for a vector of dimension dim;
 *   wrkvecimlololo is work space allocated for a vector of dimension dim;
 *   zeroQ    if true, then Q is zero and Q must be computed;
 *   noqr     flag if true, then no qr, only when not zeroQ;
 *   vrblvl   the verbose level (0 for silent).
 *
 * ON RETURN :
 *   Qrehihihi has highest doubles of the real parts of the Q in the QR;
 *   Qrelohihi has second highest doubles of the real parts of the Q in the QR;
 *   Qrehilohi has second lowest doubles of the real parts of the Q in the QR;
 *   Qrelolohi has lowest doubles of the real parts of the Q in the QR;
 *   Qrehihilo has highest doubles of the real parts of the Q in the QR;
 *   Qrelohilo has second highest doubles of the real parts of the Q in the QR;
 *   Qrehilolo has second lowest doubles of the real parts of the Q in the QR;
 *   Qrelololo has lowest doubles of the real parts of the Q in the QR;
 *   Qimhihihi has highest doubles of the imag parts of the Q in the QR;
 *   Qimlohihi has second highest doubles of the imag parts of the Q in the QR;
 *   Qimhilohi has second lowest doubles of the imag parts of the Q in the QR;
 *   Qimlolohi has lowest doubles of the imag parts of the Q in the QR;
 *   Qimhihilo has highest doubles of the imag parts of the Q in the QR;
 *   Qimlohilo has second highest doubles of the imag parts of the Q in the QR;
 *   Qimhilolo has second lowest doubles of the imag parts of the Q in the QR;
 *   Qimlololo has lowest doubles of the imag parts of the Q in the QR;
 *   Rrehihihi has highest doubles of the real parts of the R in the QR;
 *   Rrelohihi has second highest doubles of the real parts of the R in the QR;
 *   Rrehilohi has second lowest doubles of the real parts of the R in the QR;
 *   Rrelolohi has lowest doubles of the real parts of the R in the QR;
 *   Rrehihilo has highest doubles of the real parts of the R in the QR;
 *   Rrelohilo has second highest doubles of the real parts of the R in the QR;
 *   Rrehilolo has second lowest doubles of the real parts of the R in the QR;
 *   Rrelololo has lowest doubles of the real parts of the R in the QR;
 *   Rimhihihi has highest doubles of the imag parts of the R in the QR;
 *   Rimlohihi has second highest doubles of the imag parts of the R in the QR;
 *   Rimhilohi has second lowest doubles of the imag parts of the R in the QR;
 *   Rimlolohi has lowest doubles of the imag parts of the R in the QR;
 *   Rimhihilo has highest doubles of the imag parts of the R in the QR;
 *   Rimlohilo has second highest doubles of the imag parts of the R in the QR;
 *   Rimhilolo has second lowest doubles of the imag parts of the R in the QR;
 *   Rimlololo has lowest doubles of the imag parts of the R in the QR;
 *   wrkvecrehihihi is work space used to solve the linear system;
 *   wrkvecrelohihi is work space used to solve the linear system;
 *   wrkvecrehilohi is work space used to solve the linear system;
 *   wrkvecrelolohi is work space used to solve the linear system;
 *   wrkvecrehihilo is work space used to solve the linear system;
 *   wrkvecrelohilo is work space used to solve the linear system;
 *   wrkvecrehilolo is work space used to solve the linear system;
 *   wrkvecrelololo is work space used to solve the linear system;
 *   wrkvecimhihihi is work space used to solve the linear system;
 *   wrkvecimlohihi is work space used to solve the linear system;
 *   wrkvecimhilohi is work space used to solve the linear system;
 *   wrkvecimlolohi is work space used to solve the linear system;
 *   wrkvecimhihilo is work space used to solve the linear system;
 *   wrkvecimlohilo is work space used to solve the linear system;
 *   wrkvecimhilolo is work space used to solve the linear system;
 *   wrkvecimlololo is work space used to solve the linear system;
 *   solrehihihi are the highest doubles of the real parts of the head;
 *   solrelohihi are the 2nd highest doubles of the real parts of the head;
 *   solrehilohi are the 3rd highest doubles of the real parts of the head;
 *   solrehihilo are the 4th lowest doubles of the real parts of the head;
 *   solrehihilo are the 4th lowest doubles of the real parts of the head;
 *   solrelohilo are the 3rd lowest doubles of the real parts of the head;
 *   solrehilolo are the 2nd lowest doubles of the real parts of the head;
 *   solrelololo are the lowest doubles of the real parts of the head;
 *   solimhihihi are the highest doubles of the imag parts of the head;
 *   solimlohihi are the 2nd highest doubles of the imag parts of the head;
 *   solimhilohi are the 3rd highest doubles of the imag parts of the head;
 *   solimlolohi are the 4th highest doubles of the imag parts of the head;
 *   solimhihilo are the 4th lowest doubles of the imag parts of the head;
 *   solimlohilo are the 3rd lowest doubles of the imag parts of the head;
 *   solimhilolo are the 2nd lowest doubles of the imag parts of the head;
 *   solimlololo are the lowest doubles of the imag parts of the head;
 *   zeroQ    false if Q was computed;
 *   noqr      updated flag if ||dx_0|| is zero for the first time. */

void CPU_dbl8_qrbs_tail
 ( int dim, int degp1, int tailidx,
   double ***mathihihi, double ***matlohihi,
   double ***mathilohi, double ***matlolohi,
   double ***mathihilo, double ***matlohilo,
   double ***mathilolo, double ***matlololo,
   double **rhshihihi, double **rhslohihi,
   double **rhshilohi, double **rhslolohi,
   double **rhshihilo, double **rhslohilo,
   double **rhshilolo, double **rhslololo,
   double **solhihihi, double **sollohihi,
   double **solhilohi, double **sollolohi,
   double **solhihilo, double **sollohilo,
   double **solhilolo, double **sollololo,
   double **Qhihihi, double **Qlohihi, double **Qhilohi, double **Qlolohi,
   double **Qhihilo, double **Qlohilo, double **Qhilolo, double **Qlololo,
   double **Rhihihi, double **Rlohihi, double **Rhilohi, double **Rlolohi,
   double **Rhihilo, double **Rlohilo, double **Rhilolo, double **Rlololo,
   double *wrkvechihihi, double *wrkveclohihi,
   double *wrkvechilohi, double *wrkveclolohi,
   double *wrkvechihilo, double *wrkveclohilo,
   double *wrkvechilolo, double *wrkveclololo,
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
 *   mathihihi are degp1 matrices of dimension dim;
 *   matlohihi are degp1 matrices of dimension dim;
 *   mathilohi are degp1 matrices of dimension dim;
 *   matlolohi are degp1 matrices of dimension dim;
 *   mathihilo are degp1 matrices of dimension dim;
 *   matlohilo are degp1 matrices of dimension dim;
 *   mathilolo are degp1 matrices of dimension dim;
 *   matlololo are degp1 matrices of dimension dim;
 *   rhshihihi are degp1 vectors of dimension dim;
 *   rhslohihi are degp1 vectors of dimension dim;
 *   rhshilohi are degp1 vectors of dimension dim;
 *   rhslolohi are degp1 vectors of dimension dim;
 *   rhshihilo are degp1 vectors of dimension dim;
 *   rhslohilo are degp1 vectors of dimension dim;
 *   rhshilolo are degp1 vectors of dimension dim;
 *   rhslololo are degp1 vectors of dimension dim;
 *   solhihihi has space allocated for degp1 vectors of dimension dim,
 *            with the leading highest doubles defined;
 *   sollohihi  space allocated for degp1 vectors of dimension dim,
 *            with the leading second highest doubles defined;
 *   solhilohi  space allocated for degp1 vectors of dimension dim,
 *            with the leading third highest doubles defined;
 *   sollolohi  space allocated for degp1 vectors of dimension dim,
 *            with the leading fourth highest doubles defined;
 *   solhihilo  space allocated for degp1 vectors of dimension dim,
 *            with the leading fourth lowest doubles defined;
 *   sollohilo  space allocated for degp1 vectors of dimension dim,
 *            with the leading third lowest doubles defined;
 *   solhilolo  space allocated for degp1 vectors of dimension dim,
 *            with the leading second lowest doubles defined;
 *   sollololo  space allocated for degp1 vectors of dimension dim,
 *            with the leading lowest doubles defined;
 *   Qhihihi  highest doubles of the Q of the QR factorization;
 *   Qlohihi  second highest doubles of the Q of the QR factorization;
 *   Qhilohi  third highest doubles of the Q of the QR factorization;
 *   Qlolohi  fourth highest doubles of the Q of the QR factorization;
 *   Qhihilo  fourth lowest doubles of the Q of the QR factorization;
 *   Qlohilo  third lowest doubles of the Q of the QR factorization;
 *   Qhilolo  second lowest doubles of the Q of the QR factorization;
 *   Qlololo  lowest doubles of the Q of the QR factorization;
 *   Rhihihi  highest doubles of the R of the QR factorization;
 *   Rlohihi  second highest doubles of the R of the QR factorization;
 *   Rhilohi  third highest doubles of the R of the QR factorization;
 *   Rlolohi  fourth highest doubles of the R of the QR factorization;
 *   Rhihilo  fourth lowhest doubles of the R of the QR factorization;
 *   Rlohilo  third lowest doubles of the R of the QR factorization;
 *   Rhilolo  second lowest doubles of the R of the QR factorization;
 *   Rlololo  lowest doubles of the R of the QR factorization;
 *   wrkvechihihi is work space vector of dimension dim for the substitution;
 *   wrkveclohihi is work space vector of dimension dim for the substitution;
 *   wrkvechilohi is work space vector of dimension dim for the substitution;
 *   wrkveclolohi is work space vector of dimension dim for the substitution;
 *   wrkvechihilo is work space vector of dimension dim for the substitution;
 *   wrkveclohilo is work space vector of dimension dim for the substitution;
 *   wrkvechilolo is work space vector of dimension dim for the substitution;
 *   wrkveclololo is work space vector of dimension dim for the substitution;
 *   vrblvl   the verbose level (0 for silent).
 *
 * ON RETURN :
 *   rhshihihi are the updated highest doubles of right hand side;
 *   rhslohihi are the updated second highest doubles of right hand side;
 *   rhshilohi are the updated third highest doubles of right hand side;
 *   rhslolohi are the updated fourth highest doubles of right hand side;
 *   rhshihilo are the updated fourth lowest doubles of right hand side;
 *   rhslohilo are the updated third lowest doubles of right hand side;
 *   rhshilolo are the updated second lowest doubles of right hand side;
 *   rhslololo are the updated lowest doubles of right hand side;
 *   solhihihi are the highest doubles of the solution;
 *   sollohihi are the second highest doubles of the solution;
 *   solhilohi are the third highest doubles of the solution;
 *   sollolohi are the fourth highest doubles of the solution;
 *   solhihilo are the fourth lowest doubles of the solution;
 *   sollohilo are the third lowhest doubles of the solution;
 *   solhilolo are the second lowest doubles of the solution;
 *   sollololo are the lowest doubles of the solution;
 *   upidx    counts the number of updates skipped;
 *   bsidx    counts the number of backsubstitutions skipped;
 *   newtail  the new value for tailidx. */

void CPU_cmplx8_qrbs_tail
 ( int dim, int degp1, int tailidx,
   double ***matrehihihi, double ***matrelohihi,
   double ***matrehilohi, double ***matrelolohi,
   double ***matrehihilo, double ***matrelohilo,
   double ***matrehilolo, double ***matrelololo,
   double ***matimhihihi, double ***matimlohihi,
   double ***matimhilohi, double ***matimlolohi,
   double ***matimhihilo, double ***matimlohilo,
   double ***matimhilolo, double ***matimlololo,
   double **rhsrehihihi, double **rhsrelohihi,
   double **rhsrehilohi, double **rhsrelolohi,
   double **rhsrehihilo, double **rhsrelohilo,
   double **rhsrehilolo, double **rhsrelololo,
   double **rhsimhihihi, double **rhsimlohihi,
   double **rhsimhilohi, double **rhsimlolohi,
   double **rhsimhihilo, double **rhsimlohilo,
   double **rhsimhilolo, double **rhsimlololo,
   double **solrehihihi, double **solrelohihi,
   double **solrehilohi, double **solrelolohi,
   double **solrehihilo, double **solrelohilo,
   double **solrehilolo, double **solrelololo,
   double **solimhihihi, double **solimlohihi,
   double **solimhilohi, double **solimlolohi,
   double **solimhihilo, double **solimlohilo,
   double **solimhilolo, double **solimlololo,
   double **Qrehihihi, double **Qrelohihi,
   double **Qrehilohi, double **Qrelolohi,
   double **Qrehihilo, double **Qrelohilo,
   double **Qrehilolo, double **Qrelololo,
   double **Qimhihihi, double **Qimlohihi,
   double **Qimhilohi, double **Qimlolohi,
   double **Qimhihilo, double **Qimlohilo,
   double **Qimhilolo, double **Qimlololo,
   double **Rrehihihi, double **Rrelohihi,
   double **Rrehilohi, double **Rrelolohi, 
   double **Rrehihilo, double **Rrelohilo,
   double **Rrehilolo, double **Rrelololo, 
   double **Rimhihihi, double **Rimlohihi,
   double **Rimhilohi, double **Rimlolohi,
   double **Rimhihilo, double **Rimlohilo,
   double **Rimhilolo, double **Rimlololo,
   double *wrkvecrehihihi, double *wrkvecrelohihi,
   double *wrkvecrehilohi, double *wrkvecrelolohi,
   double *wrkvecrehihilo, double *wrkvecrelohilo,
   double *wrkvecrehilolo, double *wrkvecrelololo,
   double *wrkvecimhihihi, double *wrkvecimlohihi,
   double *wrkvecimhilohi, double *wrkvecimlolohi,
   double *wrkvecimhihilo, double *wrkvecimlohilo,
   double *wrkvecimhilolo, double *wrkvecimlololo,
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
 *   matrehihihi are degp1 matrices of dimension dim;
 *   matrelohihi are degp1 matrices of dimension dim;
 *   matrehilohi are degp1 matrices of dimension dim;
 *   matrelolohi are degp1 matrices of dimension dim;
 *   matrehihilo are degp1 matrices of dimension dim;
 *   matrelohilo are degp1 matrices of dimension dim;
 *   matrehilolo are degp1 matrices of dimension dim;
 *   matrelololo are degp1 matrices of dimension dim;
 *   matimhihihi are degp1 matrices of dimension dim;
 *   matimlohihi are degp1 matrices of dimension dim;
 *   matimhilohi are degp1 matrices of dimension dim;
 *   matimlolohi are degp1 matrices of dimension dim;
 *   matimhihilo are degp1 matrices of dimension dim;
 *   matimlohilo are degp1 matrices of dimension dim;
 *   matimhilolo are degp1 matrices of dimension dim;
 *   matimlololo are degp1 matrices of dimension dim;
 *   rhsrehihihi are degp1 vectors of dimension dim;
 *   rhsrelohihi are degp1 vectors of dimension dim;
 *   rhsrehilohi are degp1 vectors of dimension dim;
 *   rhsrelolohi are degp1 vectors of dimension dim;
 *   rhsrehihilo are degp1 vectors of dimension dim;
 *   rhsrelohilo are degp1 vectors of dimension dim;
 *   rhsrehilolo are degp1 vectors of dimension dim;
 *   rhsrelololo are degp1 vectors of dimension dim;
 *   rhsimhihihi are degp1 vectors of dimension dim;
 *   rhsimlohihi are degp1 vectors of dimension dim;
 *   rhsimhilohi are degp1 vectors of dimension dim;
 *   rhsimlolohi are degp1 vectors of dimension dim;
 *   rhsimhihilo are degp1 vectors of dimension dim;
 *   rhsimlohilo are degp1 vectors of dimension dim;
 *   rhsimhilolo are degp1 vectors of dimension dim;
 *   rhsimlololo are degp1 vectors of dimension dim;
 *   solrehihihi has space for degp1 vectors of dimension dim, with the
 *            highest doubles of the real parts of the head;
 *   solrelohihi has space for degp1 vectors of dimension dim, with the
 *            second highest doubles of the real parts of the head;
 *   solrehilohi has space for degp1 vectors of dimension dim, with the
 *            third highest doubles of the real parts of the head;
 *   solrelolohi has space for degp1 vectors of dimension dim, with the
 *            fourth highest doubles of the real parts of the head;
 *   solrehihilo has space for degp1 vectors of dimension dim, with the
 *            fourth lowest doubles of the real parts of the;
 *   solrelohilo has space for degp1 vectors of dimension dim, with the
 *            third lowest doubles of the real parts of the;
 *   solrehilolo has space for degp1 vectors of dimension dim, with the
 *            second lowest doubles of the real parts of the;
 *   solrelololo has space for degp1 vectors of dimension dim, with the
 *            lowest doubles of the real parts of the;
 *   solimhihihi has space for degp1 vectors of dimension dim, with the
 *            highest doubles of the imaginary parts of the head;
 *   solimlohihi has space for degp1 vectors of dimension dim, with the
 *            second highest doubles of the imaginary parts of the head;
 *   solimhilohi has space for degp1 vectors of dimension dim, with the
 *            third highest doubles of the imaginary parts of the head;
 *   solimlolohi has space for degp1 vectors of dimension dim, with the
 *            fourth highest doubles of the imaginary parts of the head;
 *   solimhihilo has space for degp1 vectors of dimension dim, with the
 *            fourth lowest doubles of the imaginary parts of the head;
 *   solimlohilo has space for degp1 vectors of dimension dim, with the
 *            third lowest doubles of the imaginary parts of the head;
 *   solimhilolo has space for degp1 vectors of dimension dim, with the
 *            second lowest doubles of the imaginary parts of the head;
 *   solimlololo has space for degp1 vectors of dimension dim, with the
 *            lowest doubles of the imaginary parts of the head;
 *   Qrehihihi has the highest doubles of the real parts of the Q of the QR;
 *   Qrelohihi has the 2nd highest doubles of the real parts of Q;
 *   Qrehilohi has the 3rd highest doubles of the real parts of Q;
 *   Qrelolohi has the 4th highest doubles of the real parts of Q;
 *   Qrehihilo has the 4th lowest doubles of the real parts of Q;
 *   Qrelohilo has the 3rd lowest doubles of the real parts of Q;
 *   Qrehilolo has the 2nd lowest doubles of the real parts of Q;
 *   Qrelololo has the lowest doubles of the real parts of Q;
 *   Qimhihihi has the highest doubles of the imaginary parts of Q;
 *   Qimlohihi has the 2nd highest doubles of the imaginary parts of Q;
 *   Qimhilohi has the 3rd highest doubles of the imaginary parts of Q;
 *   Qimlolohi has the 4th highest doubles of the imaginary parts of Q;
 *   Qimhihilo has the 4th lowest doubles of the imaginary parts of Q;
 *   Qimlohilo has the 3rd lowest doubles of the imaginary parts of Q;
 *   Qimhilolo has the 2nd lowest doubles of the imaginary parts of Q;
 *   Qimlololo has the lowest doubles of the imaginary parts of Q;
 *   Rrehihihi has the highest doubles of the real parts of the R of the QR;
 *   Rrelohihi has the 2nd highest doubles of the real parts of R;
 *   Rrehilohi has the 3rd highest doubles of the real parts of R;
 *   Rrelolohi has the 4th highest doubles of the real parts of R;
 *   Rrehihilo has the 4th lowest doubles of the real parts of R;
 *   Rrelohilo has the 3rd lowest doubles of the real parts of R;
 *   Rrehilolo has the 2nd lowest doubles of the real parts of R;
 *   Rrelololo has the lowest doubles of the real parts of R;
 *   Rimhihihi has the highest doubles of the imaginary parts of R;
 *   Rimlohihi has the 2nd highest doubles of the imaginary parts of R;
 *   Rimhilohi has the 3rd highest doubles of the imaginary parts of R;
 *   Rimlolohi has the 4th highest doubles of the imaginary parts of R;
 *   Rimhihilo has the 4th lowest doubles of the imaginary parts of R;
 *   Rimlohilo has the 3rd lowest doubles of the imaginary parts of R;
 *   Rimhilolo has the 2nd lowest doubles of the imaginary parts of R;
 *   Rimlololo has the lowest doubles of the imaginary parts of R;
 *   wrkvecrehihihi is work space of dimension dim for the substitution;
 *   wrkvecrelohihi is work space of dimension dim for the substitution;
 *   wrkvecrehilohi is work space of dimension dim for the substitution;
 *   wrkvecrelolohi is work space of dimension dim for the substitution;
 *   wrkvecrehihilo is work space of dimension dim for the substitution;
 *   wrkvecrelohilo is work space of dimension dim for the substitution;
 *   wrkvecrehilolo is work space of dimension dim for the substitution;
 *   wrkvecrelololo is work space of dimension dim for the substitution;
 *   wrkvecimhihihi is work space of dimension dim for the substitution;
 *   wrkvecimlohihi is work space of dimension dim for the substitution;
 *   wrkvecimhilohi is work space of dimension dim for the substitution;
 *   wrkvecimlolohi is work space of dimension dim for the substitution;
 *   wrkvecimhihilo is work space of dimension dim for the substitution;
 *   wrkvecimlohilo is work space of dimension dim for the substitution;
 *   wrkvecimhilolo is work space of dimension dim for the substitution;
 *   wrkvecimlololo is work space of dimension dim for the substitution;
 *   vrblvl   the verbose level (0 for silent).
 *
 * ON RETURN :
 *   rhsrehihihi are the highest doubles of the real parts
 *            of the updated right hand side;
 *   rhsrelohihi are the second highest doubles of the real parts
 *            of the updated right hand side;
 *   rhsrehilohi are the third highest doubles of the real parts
 *            of the updated right hand side;
 *   rhsrelolohi are the fourthhighest doubles of the real parts
 *            of the updated right hand side;
 *   rhsimhihilo are the fourth lowest doubles of the imaginary parts
 *            of the updated right hand side;
 *   rhsimlohilo are the third lowest doubles of the imaginary parts
 *            of the updated right hand side;
 *   rhsimhilolo are the second lowest doubles of the imaginary parts
 *            of the updated right hand side;
 *   rhsimlololo are the lowest doubles of the imaginary parts
 *            of the updated right hand side;
 *   solrehihihi are the highest doubles of the real parts of the solution;
 *   solrelohihi are the 2nd highest doubles of the real parts of sol;
 *   solrehilohi are the 3rd highest doubles of the real parts of sol;
 *   solrelolohi are the 4th highest doubles of the real parts of sol;
 *   solrehihilo are the 4th lowest doubles of the real parts of sol;
 *   solrelohilo are the 3rd lowest doubles of the real parts of sol;
 *   solrehilolo are the 2nd lowest doubles of the real parts of sol;
 *   solrelololo are the lowest doubles of the real parts of sol;
 *   solimhihihi are the highest doubles of the imaginary parts of sol;
 *   solimlohihi are the 2nd highest doubles of the imaginary parts of sol;
 *   solimhilohi are the 3rd highest doubles of the imaginary parts of sol;
 *   solimlolohi are the 4th highest doubles of the imaginary parts of sol;
 *   solimhihilo are the 4th lowest doubles of the imaginary parts of sol;
 *   solimlohilo are the 3rd lowest doubles of the imaginary parts of sol;
 *   solimhilolo are the 2nd lowest doubles of the imaginary parts of sol;
 *   solimlololo are the lowest doubles of the imaginary parts of sol;
 *   upidx    counts the number of updates skipped;
 *   bsidx    counts the number of backsubstitutions skipped;
 *   newtail  is the new value for tailidx. */

void CPU_dbl8_qrbs_solve
 ( int dim, int degp1, int tailidx,
   double ***mathihihi, double ***matlohihi,
   double ***mathilohi, double ***matlolohi,
   double ***mathihilo, double ***matlohilo,
   double ***mathilolo, double ***matlololo,
   double **rhshihihi, double **rhslohihi,
   double **rhshilohi, double **rhslolohi,
   double **rhshihilo, double **rhslohilo,
   double **rhshilolo, double **rhslololo,
   double **solhihihi, double **sollohihi,
   double **solhilohi, double **sollolohi,
   double **solhihilo, double **sollohilo,
   double **solhilolo, double **sollololo,
   double **Qhihihi, double **Qlohihi, double **Qhilohi, double **Qlolohi,
   double **Qhihilo, double **Qlohilo, double **Qhilolo, double **Qlololo,
   double **Rhihihi, double **Rlohihi, double **Rhilohi, double **Rlolohi,
   double **Rhihilo, double **Rlohilo, double **Rhilolo, double **Rlololo,
   double *wrkvechihihi, double *wrkveclohihi,
   double *wrkvechilohi, double *wrkveclolohi,
   double *wrkvechihilo, double *wrkveclohilo,
   double *wrkvechilolo, double *wrkveclololo,
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
 *   mathihihi are degp1 matrices of dimension dim;
 *   matlohihi are degp1 matrices of dimension dim;
 *   mathilohi are degp1 matrices of dimension dim;
 *   matlolohi are degp1 matrices of dimension dim;
 *   mathihilo are degp1 matrices of dimension dim;
 *   matlohilo are degp1 matrices of dimension dim;
 *   mathilolo are degp1 matrices of dimension dim;
 *   matlololo are degp1 matrices of dimension dim;
 *   rhshihihi are degp1 vectors of dimension dim;
 *   rhslohihi are degp1 vectors of dimension dim;
 *   rhshilohi are degp1 vectors of dimension dim;
 *   rhslolohi are degp1 vectors of dimension dim;
 *   rhshihilo are degp1 vectors of dimension dim;
 *   rhslohilo are degp1 vectors of dimension dim;
 *   rhshilolo are degp1 vectors of dimension dim;
 *   rhslololo are degp1 vectors of dimension dim;
 *   solhihihi has space allocated for degp1 vectors of dimension dim;
 *   sollohihi has space allocated for degp1 vectors of dimension dim;
 *   solhilohi has space allocated for degp1 vectors of dimension dim;
 *   sollolohi has space allocated for degp1 vectors of dimension dim;
 *   solhihilo has space allocated for degp1 vectors of dimension dim;
 *   sollohilo has space allocated for degp1 vectors of dimension dim;
 *   solhilolo has space allocated for degp1 vectors of dimension dim;
 *   sollololo has space allocated for degp1 vectors of dimension dim;
 *   Qhihihi  space allocated for a matrix of dimension dim;
 *   Qlohihi  space allocated for a matrix of dimension dim;
 *   Qhilohi  space allocated for a matrix of dimension dim;
 *   Qlolohi  space allocated for a matrix of dimension dim;
 *   Qhihilo  space allocated for a matrix of dimension dim;
 *   Qlohilo  space allocated for a matrix of dimension dim;
 *   Qhilolo  space allocated for a matrix of dimension dim;
 *   Qlololo  space allocated for a matrix of dimension dim;
 *   Rhihihi  space allocated for a matrix of dimension dim;
 *   Rlohihi  space allocated for a matrix of dimension dim;
 *   Rhilohi  space allocated for a matrix of dimension dim;
 *   Rlolohi  space allocated for a matrix of dimension dim;
 *   Rhihilo  space allocated for a matrix of dimension dim;
 *   Rlohilo  space allocated for a matrix of dimension dim;
 *   Rhilolo  space allocated for a matrix of dimension dim;
 *   Rlololo  space allocated for a matrix of dimension dim;
 *   wrkvechihihi has work space allocated for a vector of dimension dim;
 *   wrkveclohihi has work space allocated for a vector of dimension dim;
 *   wrkvechilohi has work space allocated for a vector of dimension dim;
 *   wrkveclolohi has work space allocated for a vector of dimension dim;
 *   wrkvechihilo has work space allocated for a vector of dimension dim;
 *   wrkveclohilo has work space allocated for a vector of dimension dim;
 *   wrkvechilolo has work space allocated for a vector of dimension dim;
 *   wrkveclololo has work space allocated for a vector of dimension dim;
 *   zeroQ    if true, then Q is zero and Q must be computed;
 *   noqr     flag if true, then no qr, only when not zeroQ;
 *   vrblvl   the verbose level (0 for silent).
 *
 * ON RETURN :
 *   Qhihihi  highest doubles of the Q of the QR of the Jacobian;
 *   Qlohihi  second highest doubles of the Q of the QR of the Jacobian;
 *   Qhilohi  third highest doubles of the Q of the QR of the Jacobian;
 *   Qlolohi  fourth highest doubles of the Q of the QR of the Jacobian;
 *   Qhihilo  fourth lowest doubles of the Q of the QR of the Jacobian;
 *   Qlohilo  third lowest doubles of the Q of the QR of the Jacobian;
 *   Qhilolo  second lowest doubles of the Q of the QR of the Jacobian;
 *   Qlololo  lowest doubles of the Q of the QR of the Jacobian;
 *   Rhihihi  highest doubles of the R in the QR of the Jacobian;
 *   Rlohihi  second highest doubles of the R in the QR of the Jacobian;
 *   Rhilohi  third highest doubles of the R in the QR of the Jacobian;
 *   Rlolohi  fourth highest doubles of the R in the QR of the Jacobian;
 *   Rhihilo  fourth lowest doubles of the R in the QR of the Jacobian;
 *   Rlohilo  third lowest doubles of the R in the QR of the Jacobian;
 *   Rhilolo  second lowest doubles of the R in the QR of the Jacobian;
 *   Rlololo  lowest doubles of the R in the QR of the Jacobian;
 *   wrkvechihihi is work space used to solve the linear system;
 *   wrkveclohihi is work space used to solve the linear system;
 *   wrkvechilohi is work space used to solve the linear system;
 *   wrkveclolohi is work space used to solve the linear system;
 *   wrkvechihilo is work space used to solve the linear system;
 *   wrkveclohilo is work space used to solve the linear system;
 *   wrkvechilolo is work space used to solve the linear system;
 *   wrkveclololo is work space used to solve the linear system;
 *   solhihihi are the highest double coefficients of the solution;
 *   sollohihi are the second highest double coefficients of the solution;
 *   solhilohi are the third highest double coefficients of the solution;
 *   sollolohi are the fourth highest double coefficients of the solution;
 *   solhihilo are the fourth lowest double coefficients of the solution;
 *   sollohilo are the third lowest double coefficients of the solution;
 *   solhilolo are the second lowest double coefficients of the solution;
 *   sollololo are the lowest double coefficients of the solution;
 *   zeroQ    false if Q was computed;
 *   noqr     flag updated when ||dx_0|| is zero for the first time;
 *   upidx    counts the number of updates skipped;
 *   bsidx    counts the number of backsubstitutions skipped;
 *   newtail  is the new value for tailidx. */

void CPU_cmplx8_qrbs_solve
 ( int dim, int degp1, int tailidx,
   double ***matrehihihi, double ***matrelohihi,
   double ***matrehilohi, double ***matrelolohi,
   double ***matrehihilo, double ***matrelohilo,
   double ***matrehilolo, double ***matrelololo,
   double ***matimhihihi, double ***matimlohihi, 
   double ***matimhilohi, double ***matimlolohi, 
   double ***matimhihilo, double ***matimlohilo, 
   double ***matimhilolo, double ***matimlololo, 
   double **rhsrehihihi, double **rhsrelohihi,
   double **rhsrehilohi, double **rhsrelolohi,
   double **rhsrehihilo, double **rhsrelohilo,
   double **rhsrehilolo, double **rhsrelololo,
   double **rhsimhihihi, double **rhsimlohihi,
   double **rhsimhilohi, double **rhsimlolohi,
   double **rhsimhihilo, double **rhsimlohilo,
   double **rhsimhilolo, double **rhsimlololo,
   double **solrehihihi, double **solrelohihi,
   double **solrehilohi, double **solrelolohi,
   double **solrehihilo, double **solrelohilo,
   double **solrehilolo, double **solrelololo,
   double **solimhihihi, double **solimlohihi,
   double **solimhilohi, double **solimlolohi,
   double **solimhihilo, double **solimlohilo,
   double **solimhilolo, double **solimlololo,
   double **Qrehihihi, double **Qrelohihi,
   double **Qrehilohi, double **Qrelolohi,
   double **Qrehihilo, double **Qrelohilo,
   double **Qrehilolo, double **Qrelololo,
   double **Qimhihihi, double **Qimlohihi,
   double **Qimhilohi, double **Qimlolohi,
   double **Qimhihilo, double **Qimlohilo,
   double **Qimhilolo, double **Qimlololo,
   double **Rrehihihi, double **Rrelohihi,
   double **Rrehilohi, double **Rrelolohi,
   double **Rrehihilo, double **Rrelohilo,
   double **Rrehilolo, double **Rrelololo,
   double **Rimhihihi, double **Rimlohihi,
   double **Rimhilohi, double **Rimlolohi,
   double **Rimhihilo, double **Rimlohilo,
   double **Rimhilolo, double **Rimlololo,
   double *wrkvecrehihihi, double *wrkvecrelohihi,
   double *wrkvecrehilohi, double *wrkvecrelolohi,
   double *wrkvecrehihilo, double *wrkvecrelohilo,
   double *wrkvecrehilolo, double *wrkvecrelololo,
   double *wrkvecimhihihi, double *wrkvecimlohihi,
   double *wrkvecimhilohi, double *wrkvecimlolohi,
   double *wrkvecimhihilo, double *wrkvecimlohilo,
   double *wrkvecimhilolo, double *wrkvecimlololo,
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
 *   matrehihihi are degp1 matrices of dimension dim;
 *   matrelohihi are degp1 matrices of dimension dim;
 *   matrehilohi are degp1 matrices of dimension dim;
 *   matrelolohi are degp1 matrices of dimension dim;
 *   matrehihilo are degp1 matrices of dimension dim;
 *   matrelohilo are degp1 matrices of dimension dim;
 *   matrehilolo are degp1 matrices of dimension dim;
 *   matrelololo are degp1 matrices of dimension dim;
 *   matimhihihi are degp1 matrices of dimension dim;
 *   matimlohihi are degp1 matrices of dimension dim;
 *   matimhilohi are degp1 matrices of dimension dim;
 *   matimlolohi are degp1 matrices of dimension dim;
 *   matimhihilo are degp1 matrices of dimension dim;
 *   matimlohilo are degp1 matrices of dimension dim;
 *   matimhilolo are degp1 matrices of dimension dim;
 *   matimlololo are degp1 matrices of dimension dim;
 *   rhsrehihihi are degp1 vectors of dimension dim;
 *   rhsrelohihi are degp1 vectors of dimension dim;
 *   rhsrehilohi are degp1 vectors of dimension dim;
 *   rhsrelolohi are degp1 vectors of dimension dim;
 *   rhsrehihilo are degp1 vectors of dimension dim;
 *   rhsrelohilo are degp1 vectors of dimension dim;
 *   rhsrehilolo are degp1 vectors of dimension dim;
 *   rhsrelololo are degp1 vectors of dimension dim;
 *   rhsimhihihi are degp1 vectors of dimension dim;
 *   rhsimlohihi are degp1 vectors of dimension dim;
 *   rhsimhilohi are degp1 vectors of dimension dim;
 *   rhsimlolohi are degp1 vectors of dimension dim;
 *   rhsimhihilo are degp1 vectors of dimension dim;
 *   rhsimlohilo are degp1 vectors of dimension dim;
 *   rhsimhilolo are degp1 vectors of dimension dim;
 *   rhsimlololo are degp1 vectors of dimension dim;
 *   solrehihihi has space allocated for degp1 vectors of dimension dim;
 *   solrelohihi has space allocated for degp1 vectors of dimension dim;
 *   solrehilohi has space allocated for degp1 vectors of dimension dim;
 *   solrelolohi has space allocated for degp1 vectors of dimension dim;
 *   solrehihilo has space allocated for degp1 vectors of dimension dim;
 *   solrelohilo has space allocated for degp1 vectors of dimension dim;
 *   solrehilolo has space allocated for degp1 vectors of dimension dim;
 *   solrelololo has space allocated for degp1 vectors of dimension dim;
 *   solimhihihi has space allocated for degp1 vectors of dimension dim;
 *   solimlohihi has space allocated for degp1 vectors of dimension dim;
 *   solimhilohi has space allocated for degp1 vectors of dimension dim;
 *   solimlolohi has space allocated for degp1 vectors of dimension dim;
 *   solimhihilo has space allocated for degp1 vectors of dimension dim;
 *   solimlohilo has space allocated for degp1 vectors of dimension dim;
 *   solimhilolo has space allocated for degp1 vectors of dimension dim;
 *   solimlololo has space allocated for degp1 vectors of dimension dim;
 *   Qrehihihi has space allocated for a matrix of dimension dim;
 *   Qrelohihi has space allocated for a matrix of dimension dim;
 *   Qrehilohi has space allocated for a matrix of dimension dim;
 *   Qrelolohi has space allocated for a matrix of dimension dim;
 *   Qrehihilo has space allocated for a matrix of dimension dim;
 *   Qrelohilo has space allocated for a matrix of dimension dim;
 *   Qrehilolo has space allocated for a matrix of dimension dim;
 *   Qrelololo has space allocated for a matrix of dimension dim;
 *   Qimhihihi has space allocated for a matrix of dimension dim;
 *   Qimlohihi has space allocated for a matrix of dimension dim;
 *   Qimhilohi has space allocated for a matrix of dimension dim;
 *   Qimlolohi has space allocated for a matrix of dimension dim;
 *   Qimhihilo has space allocated for a matrix of dimension dim;
 *   Qimlohilo has space allocated for a matrix of dimension dim;
 *   Qimhilolo has space allocated for a matrix of dimension dim;
 *   Qimlololo has space allocated for a matrix of dimension dim;
 *   Rrehihihi has space allocated for a matrix of dimension dim;
 *   Rrelohihi has space allocated for a matrix of dimension dim;
 *   Rrehilohi has space allocated for a matrix of dimension dim;
 *   Rrelolohi has space allocated for a matrix of dimension dim;
 *   Rrehihilo has space allocated for a matrix of dimension dim;
 *   Rrelohilo has space allocated for a matrix of dimension dim;
 *   Rrehilolo has space allocated for a matrix of dimension dim;
 *   Rrelololo has space allocated for a matrix of dimension dim;
 *   Rimhihihi has space allocated for a matrix of dimension dim;
 *   Rimlohihi has space allocated for a matrix of dimension dim;
 *   Rimhilohi has space allocated for a matrix of dimension dim;
 *   Rimlolohi has space allocated for a matrix of dimension dim;
 *   Rimhihilo has space allocated for a matrix of dimension dim;
 *   Rimlohilo has space allocated for a matrix of dimension dim;
 *   Rimhilolo has space allocated for a matrix of dimension dim;
 *   Rimlololo has space allocated for a matrix of dimension dim;
 *   wrkvecrehihihi has work space allocated for a vector of dimension dim;
 *   wrkvecrelohihi has work space allocated for a vector of dimension dim;
 *   wrkvecrehilohi has work space allocated for a vector of dimension dim;
 *   wrkvecrelolohi has work space allocated for a vector of dimension dim;
 *   wrkvecrehihilo has work space allocated for a vector of dimension dim;
 *   wrkvecrelohilo has work space allocated for a vector of dimension dim;
 *   wrkvecrehilolo has work space allocated for a vector of dimension dim;
 *   wrkvecrelololo has work space allocated for a vector of dimension dim;
 *   wrkvecimhihihi has work space allocated for a vector of dimension dim;
 *   wrkvecimlohihi has work space allocated for a vector of dimension dim;
 *   wrkvecimhilohi has work space allocated for a vector of dimension dim;
 *   wrkvecimlolohi has work space allocated for a vector of dimension dim;
 *   wrkvecimhihilo has work space allocated for a vector of dimension dim;
 *   wrkvecimlohilo has work space allocated for a vector of dimension dim;
 *   wrkvecimhilolo has work space allocated for a vector of dimension dim;
 *   wrkvecimlololo has work space allocated for a vector of dimension dim;
 *   zeroQ    if true, then Q is zero and Q must be computed;
 *   noqr     flag if true, then no qr, only when not zeroQ;
 *   vrblvl   the verbose level (0 for silent).
 *
 * ON RETURN :
 *   Qrehihihi has the highest doubles of the real parts of the Q in the QR;
 *   Qrelohihi has the second highest doubles of the real parts of Q;
 *   Qrehilohi has the third highest doubles of the real parts of Q;
 *   Qrelolohi has the fourth highest doubles of the real parts of Q;
 *   Qrehihilo has the fourth lowest doubles of the real parts of Q;
 *   Qrelohilo has the third lowest doubles of the real parts of Q;
 *   Qrehilolo has the second lowest doubles of the real parts of Q;
 *   Qrelololo has the lowest doubles of the real parts of Q;
 *   Qimhihihi has the highest doubles of the imaginary parts of Q;
 *   Qimlohihi has the second highest doubles of the imaginary parts of Q;
 *   Qimhilohi has the third highest doubles of the imaginary parts of Q;
 *   Qimlolohi has the fourth highest doubles of the imaginary parts of Q;
 *   Qimhihilo has the fourth lowest doubles of the imaginary parts of Q;
 *   Qimlohilo has the third lowest doubles of the imaginary parts of Q;
 *   Qimhilolo has the second lowest doubles of the imaginary parts of Q;
 *   Qimlololo has the lowest doubles of the imaginary parts of the Q;
 *   Rrehihihi has the highest doubles of the real parts of the R in the QR;
 *   Rrelohihi has the second highest doubles of the real parts of R;
 *   Rrehilohi has the third highest doubles of the real parts of R;
 *   Rrelolohi has the fourth highest doubles of the real parts of R;
 *   Rrehihilo has the fourth lowest doubles of the real parts of R;
 *   Rrelohilo has the third lowest doubles of the real parts of R;
 *   Rrehilolo has the second lowest doubles of the real parts of R;
 *   Rrelololo has the lowest doubles of the real parts of R;
 *   Rimhihihi has the highest doubles of the imaginary parts of R;
 *   Rimlohihi has the second highest doubles of the imaginary parts of R;
 *   Rimhilohi has the third highest doubles of the imaginary parts of R;
 *   Rimlolohi has the fourth highest doubles of the imaginary parts of R;
 *   Rimhihilo has the fourth lowest doubles of the imaginary parts of R;
 *   Rimlohilo has the third lowest doubles of the imaginary parts of R;
 *   Rimhilolo has the second lowest doubles of the imaginary parts of R;
 *   Rimlololo has the lowest doubles of the imaginary parts of R;
 *   wrkvecrehihihi is work space used to solve the linear systems;
 *   wrkvecrelohihi is work space used to solve the linear systems;
 *   wrkvecrehilohi is work space used to solve the linear systems;
 *   wrkvecrelolohi is work space used to solve the linear systems;
 *   wrkvecrehihilo is work space used to solve the linear systems;
 *   wrkvecrelohilo is work space used to solve the linear systems;
 *   wrkvecrehilolo is work space used to solve the linear systems;
 *   wrkvecrelololo is work space used to solve the linear systems;
 *   wrkvecimhihihi is work space used to solve the linear systems;
 *   wrkvecimlohihi is work space used to solve the linear systems;
 *   wrkvecimhilohi is work space used to solve the linear systems;
 *   wrkvecimlolohi is work space used to solve the linear systems;
 *   wrkvecimhihilo is work space used to solve the linear systems;
 *   wrkvecimlohilo is work space used to solve the linear systems;
 *   wrkvecimhilolo is work space used to solve the linear systems;
 *   wrkvecimlololo is work space used to solve the linear systems;
 *   solrehihihi are the highest doubles of the real parts of the solution;
 *   solrelohihi are the 2nd highest doubles of the real parts of sol;
 *   solrehilohi are the 3rd highest doubles of the real parts of sol;
 *   solrelolohi are the 4th highest doubles of the real parts of sol;
 *   solrehihilo are the 4th lowest doubles of the real parts of sol;
 *   solrelohilo are the 3rd lowest doubles of the real parts of sol;
 *   solrehilolo are the 2nd lowest doubles of the real parts of sol;
 *   solrelololo are the lowest doubles of the real parts of sol;
 *   solimhihihi are the highest doubles of the imag parts of sol;
 *   solimlohihi are the 2nd highest doubles of the imag parts of sol;
 *   solimhilohi are the 3rd highest doubles of the imag parts of sol;
 *   solimlolohi are the 4th highest doubles of the imag parts of sol;
 *   solimhihilo are the 4th lowest doubles of the imag parts of sol;
 *   solimlohilo are the 3rd lowest doubles of the imag parts of sol;
 *   solimhilolo are the 2nd lowest doubles of the imag parts of sol;
 *   solimlololo are the lowest doubles of the imaginary parts of sol;
 *   zeroQ    false if Q was computed;
 *   noqr     flag updated when ||dx_0|| is zero for the first time;
 *   upidx    counts the number of updates skipped;
 *   bsidx    counts the number of backsubstitutions skipped;
 *   newtail  is the new value for tailidx. */

void CPU_dbl8_linear_residue
 ( int dim, int degp1, int tailidx,
   double ***mathihihi, double ***matlohihi,
   double ***mathilohi, double ***matlolohi,
   double ***mathihilo, double ***matlohilo,
   double ***mathilolo, double ***matlololo,
   double **rhshihihi, double **rhslohihi,
   double **rhshilohi, double **rhslolohi,
   double **rhshihilo, double **rhslohilo,
   double **rhshilolo, double **rhslololo,
   double **solhihihi, double **sollohihi,
   double **solhilohi, double **sollolohi,
   double **solhihilo, double **sollohilo,
   double **solhilolo, double **sollololo,
   double **resvechihihi, double **resveclohihi,
   double **resvechilohi, double **resveclolohi,
   double **resvechihilo, double **resveclohilo,
   double **resvechilolo, double **resveclololo,
   double *resmaxhihihi, double *resmaxlohihi,
   double *resmaxhilohi, double *resmaxlolohi,
   double *resmaxhihilo, double *resmaxlohilo,
   double *resmaxhilolo, double *resmaxlololo, int vrblvl );
/*
 * DESCRIPTION :
 *   Computes the residual of the linear series system, on real data.
 *
 * ON ENTRY :
 *   dim      the dimension of the matrices and vectors;
 *   degp1    degree plus one, the size of the matrix system;
 *   tailidx  the index of the start of the update in the tail;
 *   mathihihi are degp1 matrices of dimension dim;
 *   matlohihi are degp1 matrices of dimension dim;
 *   mathilohi are degp1 matrices of dimension dim;
 *   matlolohi are degp1 matrices of dimension dim;
 *   mathihilo are degp1 matrices of dimension dim;
 *   matlohilo are degp1 matrices of dimension dim;
 *   mathilolo are degp1 matrices of dimension dim;
 *   matlololo are degp1 matrices of dimension dim;
 *   rhshihihi are degp1 right hand side vectors of dimension dim;
 *   rhslohihi are degp1 right hand side vectors of dimension dim;
 *   rhshilohi are degp1 right hand side vectors of dimension dim;
 *   rhslolohi are degp1 right hand side vectors of dimension dim;
 *   rhshihilo are degp1 right hand side vectors of dimension dim;
 *   rhslohilo are degp1 right hand side vectors of dimension dim;
 *   rhshilolo are degp1 right hand side vectors of dimension dim;
 *   rhslololo are degp1 right hand side vectors of dimension dim;
 *   solhihihi are degp1 solution vectors of dimension dim;
 *   sollohihi are degp1 solution vectors of dimension dim;
 *   solhilohi are degp1 solution vectors of dimension dim;
 *   sollolohi are degp1 solution vectors of dimension dim;
 *   solhihilo are degp1 solution vectors of dimension dim;
 *   sollohilo are degp1 solution vectors of dimension dim;
 *   solhilolo are degp1 solution vectors of dimension dim;
 *   sollololo are degp1 solution vectors of dimension dim;
 *   resvechihihi has space for the residual power series;
 *   resveclohihi has space for the residual power series;
 *   resvechilohi has space for the residual power series;
 *   resveclolohi has space for the residual power series;
 *   resvechihilo has space for the residual power series;
 *   resveclohilo has space for the residual power series;
 *   resvechilolo has space for the residual power series;
 *   resveclololo has space for the residual power series;
 *   vrblvl   is the verbose level.
 *
 * ON RETURN :
 *   resvechihihi are the highest doubles of the residual series;
 *   resveclohihi are the 2nd highest doubles of the residual series;
 *   resvechilohi are the 3rd highest doubles of the residual series;
 *   resveclolohi are the 4th highest doubles of the residual series;
 *   resvechihilo are the 4th lowest doubles of the residual series;
 *   resveclohilo are the 3rd lowest doubles of the residual series;
 *   resvechilolo are the 2nd lowest doubles of the residual series;
 *   resveclololo are the lowest doubles of the residual series;
 *   resmaxhihihi is the highest double of the maximum in resvec;
 *   resmaxlohihi is the 2nd highest double of the maximum in resvec;
 *   resmaxhilohi is the 3rd highest double of the maximum in resvec;
 *   resmaxlolohi is the 4th highest double of the maximum in resvec;
 *   resmaxhihilo is the 4th lowest double of the maximum in resvec;
 *   resmaxlohilo is the 3rd lowest double of the maximum in resvec;
 *   resmaxhilolo is the 2nd lowest double of the maximum in resvec;
 *   resmaxlololo is the lowest double of the maximum in resvec. */

void CPU_cmplx8_linear_residue
 ( int dim, int degp1, int tailidx,
   double ***matrehihihi, double ***matrelohihi,
   double ***matrehilohi, double ***matrelolohi,
   double ***matrehihilo, double ***matrelohilo,
   double ***matrehilolo, double ***matrelololo,
   double ***matimhihihi, double ***matimlohihi,
   double ***matimhilohi, double ***matimlolohi,
   double ***matimhihilo, double ***matimlohilo,
   double ***matimhilolo, double ***matimlololo,
   double **rhsrehihihi, double **rhsrelohihi,
   double **rhsrehilohi, double **rhsrelolohi,
   double **rhsrehihilo, double **rhsrelohilo,
   double **rhsrehilolo, double **rhsrelololo,
   double **rhsimhihihi, double **rhsimlohihi, 
   double **rhsimhilohi, double **rhsimlolohi, 
   double **rhsimhihilo, double **rhsimlohilo, 
   double **rhsimhilolo, double **rhsimlololo, 
   double **solrehihihi, double **solrelohihi,
   double **solrehilohi, double **solrelolohi,
   double **solrehihilo, double **solrelohilo,
   double **solrehilolo, double **solrelololo,
   double **solimhihihi, double **solimlohihi,
   double **solimhilohi, double **solimlolohi,
   double **solimhihilo, double **solimlohilo,
   double **solimhilolo, double **solimlololo,
   double **resvecrehihihi, double **resvecrelohihi,
   double **resvecrehilohi, double **resvecrelolohi,
   double **resvecrehihilo, double **resvecrelohilo,
   double **resvecrehilolo, double **resvecrelololo,
   double **resvecimhihihi, double **resvecimlohihi,
   double **resvecimhilohi, double **resvecimlolohi,
   double **resvecimhihilo, double **resvecimlohilo,
   double **resvecimhilolo, double **resvecimlololo,
   double *resmaxhihihi, double *resmaxlohihi,
   double *resmaxhilohi, double *resmaxlolohi,
   double *resmaxhihilo, double *resmaxlohilo,
   double *resmaxhilolo, double *resmaxlololo, int vrblvl );
/*
 * DESCRIPTION :
 *   Computes the residual of the linear series system, on complex data.
 *
 * ON ENTRY :
 *   dim      the dimension of the matrices and vectors;
 *   degp1    degree plus one, the size of the matrix system;
 *   tailidx  the index of the start of the update in the tail;
 *   matrehihihi are degp1 matrices of dimension dim;
 *   matrelohihi are degp1 matrices of dimension dim;
 *   matrehilohi are degp1 matrices of dimension dim;
 *   matrelolohi are degp1 matrices of dimension dim;
 *   matrehihilo are degp1 matrices of dimension dim;
 *   matrelohilo are degp1 matrices of dimension dim;
 *   matrehilolo are degp1 matrices of dimension dim;
 *   matrelololo are degp1 matrices of dimension dim;
 *   matimhihihi are degp1 matrices of dimension dim;
 *   matimlohihi are degp1 matrices of dimension dim;
 *   matimhilohi are degp1 matrices of dimension dim;
 *   matimlolohi are degp1 matrices of dimension dim;
 *   matimhihilo are degp1 matrices of dimension dim;
 *   matimlohilo are degp1 matrices of dimension dim;
 *   matimhilolo are degp1 matrices of dimension dim;
 *   matimlololo are degp1 matrices of dimension dim;
 *   rhsrehihihi are degp1 right hand side vectors of dimension dim;
 *   rhsrelohihi are degp1 right hand side vectors of dimension dim;
 *   rhsrehilohi are degp1 right hand side vectors of dimension dim;
 *   rhsrelolohi are degp1 right hand side vectors of dimension dim;
 *   rhsrehihilo are degp1 right hand side vectors of dimension dim;
 *   rhsrelohilo are degp1 right hand side vectors of dimension dim;
 *   rhsrehilolo are degp1 right hand side vectors of dimension dim;
 *   rhsrelololo are degp1 right hand side vectors of dimension dim;
 *   rhsimhihihi are degp1 right hand side vectors of dimension dim;
 *   rhsimlohihi are degp1 right hand side vectors of dimension dim;
 *   rhsimhilohi are degp1 right hand side vectors of dimension dim;
 *   rhsimlolohi are degp1 right hand side vectors of dimension dim;
 *   rhsimhihilo are degp1 right hand side vectors of dimension dim;
 *   rhsimlohilo are degp1 right hand side vectors of dimension dim;
 *   rhsimhilolo are degp1 right hand side vectors of dimension dim;
 *   rhsimlololo are degp1 right hand side vectors of dimension dim;
 *   solrehihihi are degp1 solution vectors of dimension dim;
 *   solrelohihi are degp1 solution vectors of dimension dim;
 *   solrehilohi are degp1 solution vectors of dimension dim;
 *   solrelolohi are degp1 solution vectors of dimension dim;
 *   solrehihilo are degp1 solution vectors of dimension dim;
 *   solrelohilo are degp1 solution vectors of dimension dim;
 *   solrehilolo are degp1 solution vectors of dimension dim;
 *   solrelololo are degp1 solution vectors of dimension dim;
 *   solimhihihi are degp1 solution vectors of dimension dim;
 *   solimlohihi are degp1 solution vectors of dimension dim;
 *   solimhilohi are degp1 solution vectors of dimension dim;
 *   solimlolohi are degp1 solution vectors of dimension dim;
 *   solimhihilo are degp1 solution vectors of dimension dim;
 *   solimlohilo are degp1 solution vectors of dimension dim;
 *   solimhilolo are degp1 solution vectors of dimension dim;
 *   solimlololo are degp1 solution vectors of dimension dim;
 *   resvecrehihihi has space for the residual power series;
 *   resvecrelohihi has space for the residual power series;
 *   resvecrehilohi has space for the residual power series;
 *   resvecrelolohi has space for the residual power series;
 *   resvecrehihilo has space for the residual power series;
 *   resvecrelohilo has space for the residual power series;
 *   resvecrehilolo has space for the residual power series;
 *   resvecrelololo has space for the residual power series;
 *   resvecimhihihi has space for the residual power series;
 *   resvecimlohihi has space for the residual power series;
 *   resvecimhilohi has space for the residual power series;
 *   resvecimlolohi has space for the residual power series;
 *   resvecimhihilo has space for the residual power series;
 *   resvecimlohilo has space for the residual power series;
 *   resvecimhilolo has space for the residual power series;
 *   resvecimlololo has space for the residual power series;
 *   vrblvl   is the verbose level.
 *
 * ON RETURN :
 *   resvecrehihihi are the highest doubles of the real parts
 *            of the residual series;
 *   resvecrelohihi are the second highest doubles of the real parts
 *            of the residual series;
 *   resvecrehilohi are the third highest doubles of the real parts
 *            of the residual series;
 *   resvecrelolohi are the fourth highest doubles of the real parts
 *            of the residual series;
 *   resvecrehihilo are the fourth lowest doubles of the real parts
 *            of the residual series;
 *   resvecrelohilo are the third lowest doubles of the real parts
 *            of the residual series;
 *   resvecrehilolo are the second lowest doubles of the real parts
 *            of the residual series;
 *   resvecrelololo are the lowest doubles of the real parts
 *            of the residual series;
 *   resvecimhihihi are the highest doubles of the imaginary parts
 *            of the residual series;
 *   resvecimlohihi are the second highest doubles of the imaginary parts
 *            of the residual series;
 *   resvecimhilohi are the third highest doubles of the imaginary parts
 *            of the residual series;
 *   resvecimlohihi are the fourth highest doubles of the imaginary parts
 *            of the residual series;
 *   resvecimhihilo are the fourth lowest doubles of the imaginary parts
 *            of the residual series;
 *   resvecimlohilo are the third lowest doubles of the imaginary parts
 *            of the residual series;
 *   resvecimhilolo are the second lowest doubles of the imaginary parts
 *            of the residual series;
 *   resvecimlololo are the lowest doubles of the imaginary parts
 *            of the residual series;
 *   resmaxhihihi is the highest double of the max norm of the residual;
 *   resmaxlohihi is the second highest double of the max norm of the residual;
 *   resmaxhilohi is the third highest double of the max norm of the residual;
 *   resmaxlolohi is the fourth highest double of the max norm of the residual;
 *   resmaxhihilo is the fourth lowest double of the max norm of the residual;
 *   resmaxlohilo is the third lowest double of the max norm of the residual;
 *   resmaxhilolo is the second lowest double of the max norm of the residual;
 *   resmaxlololo is the lowest double of the max norm of the residual. */

#endif
