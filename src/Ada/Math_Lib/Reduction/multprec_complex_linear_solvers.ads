with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;         
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Complex_Vectors;
with Multprec_Complex_VecVecs;
with Multprec_Complex_Matrices;

package Multprec_Complex_Linear_Solvers is

-- DESCRIPTION :
--   This package offers a few routines to solve linear systems of equations.
--   The code for lufac, lufco and lusolve is a literal translation from the
--   f77-linpack code.

  procedure Scale ( a : in out Multprec_Complex_Matrices.Matrix;
                    b : in out Multprec_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Divides the ith equation in the system a*x = b by the largest
  --   element on the ith row of a, for i in a'range(1).

  -- REQUIRED : a'range(1) = b'range(1).

  function Norm1 ( a : Multprec_Complex_Matrices.Matrix )
                 return Floating_Number;
  function Norm1 ( a : Multprec_Complex_VecVecs.VecVec )
                 return Floating_Number;

  -- DESCRIPTION :
  --   Returns the 1-norm of the matrix a.

  procedure lufac ( a : in out Multprec_Complex_Matrices.Matrix;
                    n : in integer32;
                    ipvt : out Standard_Integer_Vectors.Vector;
                    info : out integer32 );

  -- DESCRIPTION :
  --   lufac factors a complex matrix by gaussian elimination

  --   lufac is usually called by lufco, but it can be called
  --   directly with a saving of time if rcond is not needed.
  --   (time for lufco) = (1 + 9/n)*(time for lufac).

  -- ON ENTRY :
  --   a       complex matrix(1..n,1..n) to be factored
  --   n       the dimension of the matrix a

  -- ON RETURN :
  --   a       an upper triangular matrix and the multipliers
  --           which were used to obtain it.
  --           The factorization can be written a = l*u where
  --           l is a product of permutation and unit lower
  --           triangular matrices and u is upper triangular.
  --   ipvt    an integer vector of pivot indices
  --   info    = 0 normal value,
  --           = k if u(k,k) = 0.0.
  --           This is not an error for this routine, but it indicates
  --           that lusolve will divide by zero if called.  Use rcond in
  --           lufco for a reliable indication of singularity.

  procedure lufac ( a : in out Multprec_Complex_VecVecs.VecVec;
                    n : in integer32;
                    ipvt : out Standard_Integer_Vectors.Vector;
                    info : out integer32 );

  -- DESCRIPTION :
  --   LU factorization on vector of vectors data type.
  --   Except for a, the parameters n, ipvt, and info play the same role
  --   as the lufac on a matrix.
  --   The columns of the matrix a are stored as vectors
  --   and the ranges of the vectors are supposed to contain 1..n.

  procedure estco ( a : in Multprec_Complex_Matrices.Matrix;
                    n : in integer32;
                    ipvt : in Standard_Integer_Vectors.Vector;
                    anorm : in Floating_Number; rcond : out Floating_Number );

  -- DESCRIPTION :
  --   estco estimates the condition number of the matrix,
  --   based on the given lu factorization of the matrix.

  -- ON ENTRY :
  --   a       an upper triangular matrix and the multipliers 
  --           which are used to obtain it.
  --           The factorization can be written a = l*u, where
  --           l is a product of permutation and unit lower triangular
  --           matrices and u is upper triangular.
  --   n       the dimension of the matrix a.
  --   ipvt    an integer vector of pivot indices.
  --   anorm   1-norm of the matrix a, computed before factorization of a.

  -- ON RETURN :
  --   rcond   an estimate of the reciprocal condition of a.
  --           For the system a*x = b, relative perturbations
  --           in a and b of size epsilon may cause relative
  --           perturbations in x of size epsilon/rcond.
  --           If rcond is so small that the logical expression
  --                  1.0 + rcond = 1.0 
  --           is true, than a may be singular to working precision.
  --           In particular, rcond is zero if exact singularity is
  --           detected or the estimate underflows.

  procedure estco ( a : in Multprec_Complex_VecVecs.VecVec;
                    n : in integer32;
                    ipvt : in Standard_Integer_Vectors.Vector;
                    anorm : in Floating_Number; rcond : out Floating_Number );

  -- DESCRIPTION :
  --   Version of estco for matrices as vectors of columns.

  procedure lufco ( a : in out Multprec_Complex_Matrices.Matrix;
                    n : in integer32;
                    ipvt : out Standard_Integer_Vectors.Vector;
                    rcond : out Floating_Number );

  -- DESCRIPTION :
  --   lufco factors a complex matrix by gaussian elimination
  --   and estimates the condition of the matrix.

  -- If rcond is not needed, lufac is slightly faster.
  -- To solve a*x = b, follow lufco by lusolve.

  -- ON ENTRY :
  --   a       complex matrix(1..n,1..n) to be factored
  --   n       the dimension of the matrix a

  -- ON RETURN :
  --   a       an upper triangular matrix and the multipliers 
  --           which are used to obtain it.
  --           The factorization can be written a = l*u, where
  --           l is a product of permutation and unit lower triangular
  --           matrices and u is upper triangular.
  --   ipvt    an integer vector of pivot indices
  --   rcond   an estimate of the reciprocal condition of a.
  --           For the system a*x = b, relative perturbations
  --           in a and b of size epsilon may cause relative
  --           perturbations in x of size epsilon/rcond.
  --           If rcond is so small that the logical expression
  --                  1.0 + rcond = 1.0 
  --           is true, than a may be singular to working precision.
  --           In particular, rcond is zero if exact singularity is
  --           detected or the estimate underflows.

  procedure lufco ( a : in out Multprec_Complex_VecVecs.VecVec;
                    n : in integer32;
                    ipvt : out Standard_Integer_Vectors.Vector;
                    rcond : out Floating_Number );

  -- DESCRIPTION :
  --   Version of lufco for matrices as vectors of columns.

  procedure lusolve ( a : in Multprec_Complex_Matrices.Matrix;
                      n : in integer32;
                      ipvt : in Standard_Integer_Vectors.Vector;
                      b : in out Multprec_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   lusolve solves the complex system a*x = b using the factors
  --   computed by lufac or lufco.

  -- ON ENTRY :
  --   a       a complex matrix(1..n,1..n), the output from lufac or lufco;
  --   n       the dimension of the matrix a;
  --   ipvt    the pivot vector from lufac or lufco;
  --   b       the right hand side vector.

  -- ON RETURN :
  --   b       the solution vector x.

  procedure lusolve ( a : in Multprec_Complex_VecVecs.VecVec;
                      n : in integer32;
                      ipvt : in Standard_Integer_Vectors.Vector;
                      b : in out Multprec_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Version of lusolve for matrices as vectors of columns.

  procedure Triangulate ( a : in out Multprec_Complex_Matrices.Matrix;
                          tol : in double_float;
                          size : in natural32; n,m : in integer32 );
  procedure Triangulate ( file : in file_type;
                          a : in out Multprec_Complex_Matrices.Matrix;
                          tol : in double_float;
                          size : in natural32; n,m : in integer32 );

  -- DESCRIPTION :
  --   Triangulate makes the n*m complex matrix a triangular using
  --   Gaussian elimination.
  --   When "file" is given as input parameter, "shadowed" calculations
  --   on standard complex numbers are performed with comparisons.

  -- ON ENTRY :
  --   file    for intermediate diagnostics;
  --   a       a complex matrix(1..n,1..m);
  --   tol     tolerance to decide whether a number is zero or not;
  --   size    size of the numbers;
  --   n       the number of rows of a;
  --   m       the number of columns of a.

  -- ON RETURN :
  --   a       the triangulated matrix.

  procedure Diagonalize ( a : in out Multprec_Complex_Matrices.Matrix;
                          n,m : in integer32 );

  -- DESCRIPTION :
  --   diagonalize makes the n*m complex matrix a diagonal using
  --   the Gauss-Jordan method.

  -- ON ENTRY :
  --   a       a complex matrix(1..n,1..m);
  --   n       the number of rows of a;
  --   m       the number of columns of a.

  -- ON RETURN :
  --   a       the diagonalized matrix.

end Multprec_Complex_Linear_Solvers;
