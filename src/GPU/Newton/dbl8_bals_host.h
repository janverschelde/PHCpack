// The file dbl8_bals_host.h specifies functions to solve linear systems
// of power series by linearization in octo double precision.

#ifndef __dbl8_bals_host_h__
#define __dbl8_bals_host_h__

void CPU_dbl8_linear_head
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

void CPU_dbl8_linear_tail
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

void CPU_dbl8_linear_solve
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
