// The file dbl4_bals_host.h specifies functions to solve linear systems
// of power series by linearization in quad double precision.

#ifndef __dbl4_bals_host_h__
#define __dbl4_bals_host_h__

void CPU_dbl4_linear_head
 ( int dim, int degp1,
   double ***mathihi, double ***matlohi, double ***mathilo, double ***matlolo, 
   double **rhshihi, double **rhslohi, double **rhshilo, double **rhslolo,
   double **solhihi, double **sollohi, double **solhilo, double **sollolo,
   double **wrkmathihi, double **wrkmatlohi,
   double **wrkmathilo, double **wrkmatlolo,
   double *wrkvechihi, double *wrkveclohi,
   double *wrkvechilo, double *wrkveclolo, int *pivots, int vrblvl );
/*
 * DESCRIPTION :
 *   Computes the leading terms of the power series solution
 *   to a linear system of power series, in linearized format.
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
 *   wrkmathihi is work space allocated for a matrix of dimension dim;
 *   wrkmatlohi is work space allocated for a matrix of dimension dim;
 *   wrkmathilo is work space allocated for a matrix of dimension dim;
 *   wrkmatlolo is work space allocated for a matrix of dimension dim;
 *   wrkvechihi is work space allocated for a vector of dimension dim;
 *   wrkveclohi is work space allocated for a vector of dimension dim;
 *   wrkvechilo is work space allocated for a vector of dimension dim;
 *   wrkveclolo is work space allocated for a vector of dimension dim;
 *   pivots   space for dim integers to store the pivots;
 *   vrblvl   the verbose level (0 for silent).
 *
 * ON RETURN :
 *   wrkmathihi are the highest doubles of the LU factorization
 *            of the Jacobian matrix;
 *   wrkmatlohi are the second highest doubles of the LU factorization
 *            of the Jacobian matrix;
 *   wrkmathilo are the second lowest doubles of the LU factorization
 *            of the Jacobian matrix;
 *   wrkmatlolo are the lowest doubles of the LU factorization
 *            of the Jacobian matrix;
 *   pivots   pivots used in the LU factorization of the Jacobian;
 *   solhihi  the highest doubles of the solution series;
 *   sollohi  the second highest doubles of the solution series;
 *   solhilo  the second lowest doubles of the solution series;
 *   sollolo  the lowest doubles of the solution series. */

void CPU_dbl4_linear_tail
 ( int dim, int degp1,
   double ***mathihi, double ***matlohi, double ***mathilo, double ***matlolo, 
   double **rhshihi, double **rhslohi, double **rhshilo, double **rhslolo,
   double **solhihi, double **sollohi, double **solhilo, double **sollolo,
   double **wrkmathihi, double **wrkmatlohi,
   double **wrkmathilo, double **wrkmatlolo, int *pivots, int vrblvl );
/*
 * DESCRIPTION :
 *   Computes the trailing terms of the power series solution
 *   to a linear system of power series, in linearized format.
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
 *   solhihi  space allocated for degp1 vectors of dimension dim,
 *            with the highest doubles of the leading coefficients defined;
 *   sollohi  space allocated for degp1 vectors of dimension dim,
 *            with the 2nd highest doubles of the leading coefficients defined;
 *   solhilo  space allocated for degp1 vectors of dimension dim,
 *            with the 2nd lowest doubles of the leading coefficients defined;
 *   sollolo  space allocated for degp1 vectors of dimension dim,
 *            with the lowest doubles of the leading coefficients defined;
 *   wrkmathihi are the highest doubles of the LU factorization
 *            of the Jacobian matrix;
 *   wrkmatlohi are the second highest doubles of the LU factorization
 *            of the Jacobian matrix;
 *   wrkmathilo are the second lowest doubles of the LU factorization
 *            of the Jacobian matrix;
 *   wrkmatlolo are the lowest doubles of the LU factorization
 *            of the Jacobian matrix;
 *   pivots   dim integers store the pivots of the LU factorization;
 *   vrblvl   the verbose level (0 for silent).
 *
 * ON RETURN :
 *   rhshihi  updated right hand side used as work space;
 *   rhslohi  updated right hand side used as work space;
 *   rhshilo  updated right hand side used as work space;
 *   rhslolo  updated right hand side used as work space;
 *   solhihi  the highest doubles of the solution series;
 *   sollohi  the second highest doubles of the solution series;
 *   solhilo  the second lowest doubles of the solution series;
 *   sollolo  the lowest doubles of the solution series. */

void CPU_dbl4_linear_solve
 ( int dim, int degp1,
   double ***mathihi, double ***matlohi, double ***mathilo, double ***matlolo, 
   double **rhshihi, double **rhslohi, double **rhshilo, double **rhslolo,
   double **solhihi, double **sollohi, double **solhilo, double **sollolo,
   double **wrkmathihi, double **wrkmatlohi,
   double **wrkmathilo, double **wrkmatlolo,
   double *wrkvechihi, double *wrkveclohi,
   double *wrkvechilo, double *wrkveclolo, int *pivots, int vrblvl );
/*
 * DESCRIPTION :
 *   Solves a linear system of power series, in linearized format.
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
 *   wrkmathihi has work space allocated for a matrix of dimension dim;
 *   wrkmatlohi has work space allocated for a matrix of dimension dim;
 *   wrkmathilo has work space allocated for a matrix of dimension dim;
 *   wrkmatlolo has work space allocated for a matrix of dimension dim;
 *   wrkvechihi has work space allocated for a vector of dimension dim;
 *   wrkveclohi has work space allocated for a vector of dimension dim;
 *   wrkvechilo has work space allocated for a vector of dimension dim;
 *   wrkveclolo has work space allocated for a vector of dimension dim;
 *   pivots   space for dim integers to store the pivots;
 *   vrblvl   the verbose level (0 for silent).
 *
 * ON RETURN :
 *   wrkmathihi are the highest doubles of the LU factorization
 *            of the Jacobian matrix;
 *   wrkmatlohi are the second highest doubles of the LU factorization
 *            of the Jacobian matrix;
 *   wrkmathilo are the second lowest doubles of the LU factorization
 *            of the Jacobian matrix;
 *   wrkmatlolo are the lowest doubles of the LU factorization
 *            of the Jacobian matrix;
 *   pivots   pivots used in the LU factorization of the Jacobian;
 *   solhihi  the highest doubles of the solution series;
 *   sollohi  the second highest doubles of the solution series;
 *   solhilo  the second lowest doubles of the solution series;
 *   sollolo  the lowest doubles of the solution series. */

void CPU_dbl4_linear_residue
 ( int dim, int degp1,
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

#endif
