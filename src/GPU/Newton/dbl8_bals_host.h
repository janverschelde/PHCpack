// The file dbl8_bals_host.h specifies functions to solve linear systems
// of power series by linearization in octo double precision.

#ifndef __dbl8_bals_host_h__
#define __dbl8_bals_host_h__

void CPU_dbl8_lusb_head
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
   double **wrkmathihihi, double **wrkmatlohihi,
   double **wrkmathilohi, double **wrkmatlolohi,
   double **wrkmathihilo, double **wrkmatlohilo,
   double **wrkmathilolo, double **wrkmatlololo,
   double *wrkvechihihi, double *wrkveclohihi,
   double *wrkvechilohi, double *wrkveclolohi,
   double *wrkvechihilo, double *wrkveclohilo,
   double *wrkvechilolo, double *wrkveclololo, int *pivots, int vrblvl );
/*
 * DESCRIPTION :
 *   Computes the leading terms of the power series solution
 *   to a linear system of power series, in linearized format.
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
 *   wrkmathihihi is work space allocated for a matrix of dimension dim;
 *   wrkmatlohihi is work space allocated for a matrix of dimension dim;
 *   wrkmathilohi is work space allocated for a matrix of dimension dim;
 *   wrkmatlolohi is work space allocated for a matrix of dimension dim;
 *   wrkmathihilo is work space allocated for a matrix of dimension dim;
 *   wrkmatlohilo is work space allocated for a matrix of dimension dim;
 *   wrkmathilolo is work space allocated for a matrix of dimension dim;
 *   wrkmatlololo is work space allocated for a matrix of dimension dim;
 *   wrkvechihihi is work space allocated for a vector of dimension dim;
 *   wrkveclohihi is work space allocated for a vector of dimension dim;
 *   wrkvechilohi is work space allocated for a vector of dimension dim;
 *   wrkveclolohi is work space allocated for a vector of dimension dim;
 *   wrkvechihilo is work space allocated for a vector of dimension dim;
 *   wrkveclohilo is work space allocated for a vector of dimension dim;
 *   wrkvechilolo is work space allocated for a vector of dimension dim;
 *   wrkveclololo is work space allocated for a vector of dimension dim;
 *   pivots   space for dim integers to store the pivots;
 *   vrblvl   the verbose level (0 for silent).
 *
 * ON RETURN :
 *   wrkmathihihi are the highest doubles of the LU factorization
 *            of the Jacobian matrix;
 *   wrkmatlohihi are the second highest doubles of the LU factorization
 *            of the Jacobian matrix;
 *   wrkmathilohi are the third highest doubles of the LU factorization
 *            of the Jacobian matrix;
 *   wrkmatlolohi are the fourth highest doubles of the LU factorization
 *            of the Jacobian matrix;
 *   wrkmathihilo are the fourth lowst doubles of the LU factorization
 *            of the Jacobian matrix;
 *   wrkmatlohilo are the third lowest doubles of the LU factorization
 *            of the Jacobian matrix;
 *   wrkmathilolo are the second lowes doubles of the LU factorization
 *            of the Jacobian matrix;
 *   wrkmatlololo are the lowest doubles of the LU factorization
 *            of the Jacobian matrix;
 *   pivots   pivots used in the LU factorization of the Jacobian;
 *   solhihihi are the highest doubles of the solution series;
 *   sollohihi are the second highest doubles of the solution series;
 *   solhilohi are the third highest doubles of the solution series;
 *   sollolohi are the fourth highest doubles of the solution series.
 *   solhihilo are the fourth lowest doubles of the solution series;
 *   sollohilo are the third lowest doubles of the solution series;
 *   solhilolo are the second lowest doubles of the solution series;
 *   sollololo are the lowest doubles of the solution series. */

void CPU_dbl8_qrbs_tail
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
   double *wrkvechilolo, double *wrkveclololo, int vrblvl );
/*
 * DESCRIPTION :
 *   Computes the trailing terms of the power series solution
 *   to a linear system of power series, in linearized format,
 *   applying substitution given a QR factorization.
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
 *   sollololo are the lowest doubles of the solution. */

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
   double **wrkmathihihi, double **wrkmatlohihi,
   double **wrkmathilohi, double **wrkmatlolohi,
   double **wrkmathihilo, double **wrkmatlohilo,
   double **wrkmathilolo, double **wrkmatlololo,
   double **Qhihihi, double **Qlohihi, double **Qhilohi, double **Qlolohi,
   double **Qhihilo, double **Qlohilo, double **Qhilolo, double **Qlololo,
   double **Rhihihi, double **Rlohihi, double **Rhilohi, double **Rlolohi, 
   double **Rhihilo, double **Rlohilo, double **Rhilolo, double **Rlololo, 
   double *wrkvechihihi, double *wrkveclohihi,
   double *wrkvechilohi, double *wrkveclolohi,
   double *wrkvechihilo, double *wrkveclohilo,
   double *wrkvechilolo, double *wrkveclololo, int vrblvl );
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
 *   wrkmathihihi is work space allocated for a matrix of dimension dim;
 *   wrkmatlohihi is work space allocated for a matrix of dimension dim;
 *   wrkmathilohi is work space allocated for a matrix of dimension dim;
 *   wrkmatlolohi is work space allocated for a matrix of dimension dim;
 *   wrkmathihilo is work space allocated for a matrix of dimension dim;
 *   wrkmatlohilo is work space allocated for a matrix of dimension dim;
 *   wrkmathilolo is work space allocated for a matrix of dimension dim;
 *   wrkmatlololo is work space allocated for a matrix of dimension dim;
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
 *   vrblvl   the verbose level (0 for silent).
 *
 * ON RETURN :
 *   wrkmathihihi has a copy of the Jacobian matrix in mathihihi[0];
 *   wrkmatlohihi has a copy of the Jacobian matrix in matlohihi[0];
 *   wrkmathilohi has a copy of the Jacobian matrix in mathilohi[0];
 *   wrkmatlolohi has a copy of the Jacobian matrix in matlolohi[0];
 *   wrkmathihilo has a copy of the Jacobian matrix in mathihilo[0];
 *   wrkmatlohilo has a copy of the Jacobian matrix in matlohilo[0];
 *   wrkmathilolo has a copy of the Jacobian matrix in mathilolo[0];
 *   wrkmatlololo has a copy of the Jacobian matrix in matlololo[0];
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
 *   sollololo  lowest doubles of the coefficients of the solution. */

void CPU_dbl8_lusb_tail
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
   double **wrkmathihihi, double **wrkmatlohihi,
   double **wrkmathilohi, double **wrkmatlolohi,
   double **wrkmathihilo, double **wrkmatlohilo,
   double **wrkmathilolo, double **wrkmatlololo, int *pivots, int vrblvl );
/*
 * DESCRIPTION :
 *   Computes the trailing terms of the power series solution
 *   to a linear system of power series, in linearized format.
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
 *   solhihihi has space allocated for degp1 vectors of dimension dim,
 *            with the highest doubles of the leading coefficients;
 *   sollohihi has space allocated for degp1 vectors of dimension dim,
 *            with the 2nd highest doubles of the leading coefficients;
 *   solhilohi has space allocated for degp1 vectors of dimension dim,
 *            with the 3rd highest doubles of the leading coefficients;
 *   sollolohi has space allocated for degp1 vectors of dimension dim,
 *            with the 4th highest doubles of the leading coefficients;
 *   solhihilo has space allocated for degp1 vectors of dimension dim,
 *            with the 4th lowest doubles of the leading coefficients;
 *   sollohilo has space allocated for degp1 vectors of dimension dim,
 *            with the 3rd lowest doubles of the leading coefficients;
 *   solhilolo has space allocated for degp1 vectors of dimension dim,
 *            with the 2nd lowest doubles of the leading coefficients;
 *   sollololo has space allocated for degp1 vectors of dimension dim,
 *            with the lowest doubles of the leading coefficients;
 *   wrkmathihihi are the highest doubles of the LU factorization
 *            of the Jacobian matrix;
 *   wrkmatlohihi are the second highest doubles of the LU factorization
 *            of the Jacobian matrix;
 *   wrkmathilohi are the third highest doubles of the LU factorization
 *            of the Jacobian matrix;
 *   wrkmatlolohi are the fourth highest doubles of the LU factorization
 *            of the Jacobian matrix;
 *   wrkmathihilo are the fourth lowest doubles of the LU factorization
 *            of the Jacobian matrix;
 *   wrkmatlohilo are the third lowest doubles of the LU factorization
 *            of the Jacobian matrix;
 *   wrkmathilolo are the second lowest doubles of the LU factorization
 *            of the Jacobian matrix;
 *   wrkmatlololo are the lowest doubles of the LU factorization
 *            of the Jacobian matrix;
 *   pivots   dim integers store the pivots of the LU factorization;
 *   vrblvl   the verbose level (0 for silent).
 *
 * ON RETURN :
 *   rhshihihi is updated right hand side used as work space;
 *   rhslohihi is updated right hand side used as work space;
 *   rhshilohi is updated right hand side used as work space;
 *   rhslolohi is updated right hand side used as work space;
 *   rhshihilo is updated right hand side used as work space;
 *   rhslohilo is updated right hand side used as work space;
 *   rhshilolo is updated right hand side used as work space;
 *   rhslololo is updated right hand side used as work space;
 *   solhihihi are the highest doubles of the solution series;
 *   sollohihi are the second highest doubles of the solution series;
 *   solhilohi are the third highest doubles of the solution series;
 *   sollolohi are the fourth highest doubles of the solution series;
 *   solhihilo are the fourth lowest doubles of the solution series;
 *   sollohilo are the third lowest doubles of the solution series;
 *   solhilolo are the second lowest doubles of the solution series;
 *   sollololo are the lowest doubles of the solution series. */

void CPU_dbl8_lusb_solve
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
   double **wrkmathihihi, double **wrkmatlohihi,
   double **wrkmathilohi, double **wrkmatlolohi,
   double **wrkmathihilo, double **wrkmatlohilo,
   double **wrkmathilolo, double **wrkmatlololo,
   double *wrkvechihihi, double *wrkveclohihi,
   double *wrkvechilohi, double *wrkveclolohi,
   double *wrkvechihilo, double *wrkveclohilo,
   double *wrkvechilolo, double *wrkveclololo, int *pivots, int vrblvl );
/*
 * DESCRIPTION :
 *   Solves a linear system of power series, in linearized format.
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
 *   solhihihi has space allocated for degp1 vectors of dimension dim;
 *   sollohihi has space allocated for degp1 vectors of dimension dim;
 *   solhilolo has space allocated for degp1 vectors of dimension dim;
 *   sollololo has space allocated for degp1 vectors of dimension dim;
 *   solhilolo has space allocated for degp1 vectors of dimension dim;
 *   sollololo has space allocated for degp1 vectors of dimension dim;
 *   wrkmathihihi has work space allocated for a matrix of dimension dim;
 *   wrkmatlohihi has work space allocated for a matrix of dimension dim;
 *   wrkmathilohi has work space allocated for a matrix of dimension dim;
 *   wrkmatlolohi has work space allocated for a matrix of dimension dim;
 *   wrkmathihilo has work space allocated for a matrix of dimension dim;
 *   wrkmatlohilo has work space allocated for a matrix of dimension dim;
 *   wrkmathilolo has work space allocated for a matrix of dimension dim;
 *   wrkmatlololo has work space allocated for a matrix of dimension dim;
 *   wrkvechihihi has work space allocated for a vector of dimension dim;
 *   wrkveclohihi has work space allocated for a vector of dimension dim;
 *   wrkvechilohi has work space allocated for a vector of dimension dim;
 *   wrkveclolohi has work space allocated for a vector of dimension dim;
 *   wrkvechihilo has work space allocated for a vector of dimension dim;
 *   wrkveclohilo has work space allocated for a vector of dimension dim;
 *   wrkvechilolo has work space allocated for a vector of dimension dim;
 *   wrkveclololo has work space allocated for a vector of dimension dim;
 *   pivots   space for dim integers to store the pivots;
 *   vrblvl   the verbose level (0 for silent).
 *
 * ON RETURN :
 *   wrkmathihihi are the highest doubles of the LU factorization
 *            of the Jacobian matrix;
 *   wrkmatlohihi are the second highest doubles of the LU factorization
 *            of the Jacobian matrix;
 *   wrkmathilohi are the third highest doubles of the LU factorization
 *            of the Jacobian matrix;
 *   wrkmatlolohi are the fourth highest doubles of the LU factorization
 *            of the Jacobian matrix;
 *   wrkmathihilo are the fourth lowest doubles of the LU factorization
 *            of the Jacobian matrix;
 *   wrkmatlohilo are the third lowest doubles of the LU factorization
 *            of the Jacobian matrix;
 *   wrkmathilolo are the second lowest doubles of the LU factorization
 *            of the Jacobian matrix;
 *   wrkmatlololo are the lowest doubles of the LU factorization
 *            of the Jacobian matrix;
 *   pivots   pivots used in the LU factorization of the Jacobian;
 *   solhihihi are the highest doubles of the solution series;
 *   sollohihi are the second highest doubles of the solution series;
 *   solhilohi are the third highest doubles of the solution series;
 *   sollolohi are the fourth highest doubles of the solution series;
 *   solhihilo are the fourth lowest doubles of the solution series;
 *   sollohilo are the third lowest doubles of the solution series;
 *   solhilolo are the second lowest doubles of the solution series;
 *   sollololo are the lowest doubles of the solution series. */

void CPU_dbl8_qrbs_solve
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
   double **wrkmathihihi, double **wrkmatlohihi,
   double **wrkmathilohi, double **wrkmatlolohi,
   double **wrkmathihilo, double **wrkmatlohilo,
   double **wrkmathilolo, double **wrkmatlololo,
   double **Qhihihi, double **Qlohihi, double **Qhilohi, double **Qlolohi,
   double **Qhihilo, double **Qlohilo, double **Qhilolo, double **Qlololo,
   double **Rhihihi, double **Rlohihi, double **Rhilohi, double **Rlolohi,
   double **Rhihilo, double **Rlohilo, double **Rhilolo, double **Rlololo,
   double *wrkvechihihi, double *wrkveclohihi,
   double *wrkvechilohi, double *wrkveclolohi,
   double *wrkvechihilo, double *wrkveclohilo,
   double *wrkvechilolo, double *wrkveclololo, int vrblvl );
/*
 * DESCRIPTION :
 *   Solves a linear system of power series, in linearized format,
 *   using QR factorization and substitutions.
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
 *   wrkmathihihi has work space allocated for a matrix of dimension dim;
 *   wrkmatlohihi has work space allocated for a matrix of dimension dim;
 *   wrkmathihihi has work space allocated for a matrix of dimension dim;
 *   wrkmatlohihi has work space allocated for a matrix of dimension dim;
 *   wrkmathilolo has work space allocated for a matrix of dimension dim;
 *   wrkmatlololo has work space allocated for a matrix of dimension dim;
 *   wrkmathilolo has work space allocated for a matrix of dimension dim;
 *   wrkmatlololo has work space allocated for a matrix of dimension dim;
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
 *   vrblvl   the verbose level (0 for silent).
 *
 * ON RETURN :
 *   wrkmathihihi has a copy of the Jacobian matrix in mathihihi[0];
 *   wrkmatlohihi has a copy of the Jacobian matrix in matlohihi[0];
 *   wrkmathilohi has a copy of the Jacobian matrix in mathilohi[0];
 *   wrkmatlolohi has a copy of the Jacobian matrix in matlolohi[0];
 *   wrkmathihilo has a copy of the Jacobian matrix in mathihilo[0];
 *   wrkmatlohilo has a copy of the Jacobian matrix in matlohilo[0];
 *   wrkmathilolo has a copy of the Jacobian matrix in mathilolo[0];
 *   wrkmatlololo has a copy of the Jacobian matrix in matlololo[0];
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
 *   sollololo are the lowest double coefficients of the solution. */

void CPU_dbl8_linear_residue
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
 *   Computes the residual of the linear power series system.
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

#endif
