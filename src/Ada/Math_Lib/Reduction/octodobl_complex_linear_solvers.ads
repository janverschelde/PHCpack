with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;         
with Standard_Natural_Matrices;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Octo_Double_Numbers;                use Octo_Double_Numbers;
with OctoDobl_Complex_Vectors;
with OctoDobl_Complex_VecVecs;
with OctoDobl_Complex_Matrices;

package OctoDobl_Complex_Linear_Solvers is

-- DESCRIPTION :
--   This package offers a few routines to solve linear systems of equations,
--   with complex arithmetic in octo double precision.
--   Two different definitions of matrices are supported:
--   (1) with the builtin two dimensional array type,
--   (2) as vectors of columns, for more efficient memory movements.
--   The code for lufac, lufco and lusolve is a literal translation from the
--   f77-linpack code.

  procedure Scale ( a : in out OctoDobl_Complex_Matrices.Matrix;
                    b : in out OctoDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Divides the ith equation in the system a*x = b by the largest
  --   element on the ith row of a, for i in a'range(1).

  -- REQUIRED : a'range(1) = b'range(1).

  function Norm1 ( a : OctoDobl_Complex_Matrices.Matrix ) return octo_double;
  function Norm1 ( a : OctoDobl_Complex_VecVecs.VecVec ) return octo_double;

  -- DESCRIPTION :
  --   Returns the 1-norm of the matrix a.

  procedure lufac ( a : in out OctoDobl_Complex_Matrices.Matrix;
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
  --   info    = 0  normal value
  --           = k  if u(k,k) = 0.0.
  --                This is not an error for this routine,
  --                but it does indicate that lusolve will
  --                divide by zero if called.  Use rcond in
  --                lufco for a reliable indication of singularity.

  procedure lufac ( a : in out OctoDobl_Complex_VecVecs.VecVec;
                    n : in integer32;
                    ipvt : out Standard_Integer_Vectors.Vector;
                    info : out integer32 );

  -- DESCRIPTION :
  --   LU factorization on vector of vectors data type.
  --   Except for a, the parameters n, ipvt, and info play the same role
  --   as the lufac on a matrix.
  --   The columns of the matrix a are stored as vectors
  --   and the ranges of the vectors are supposed to contain 1..n.

  procedure estco ( a : in OctoDobl_Complex_Matrices.Matrix;
                    n : in integer32;
                    ipvt : in Standard_Integer_Vectors.Vector;
                    anorm : in octo_double; rcond : out octo_double );

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

  procedure estco ( a : in OctoDobl_Complex_VecVecs.VecVec;
                    n : in integer32;
                    ipvt : in Standard_Integer_Vectors.Vector;
                    anorm : in octo_double; rcond : out octo_double );

  -- DESCRIPTION :
  --   Estimation of the condition number with in a the output of lufac.
  --   All parameters play the same role as the estco on a matrix.
  --   The columns of the matrix a are stored as vectors
  --   and the ranges of the vectors are supposed to contain 1..n.

  procedure lufco ( a : in out OctoDobl_Complex_Matrices.Matrix;
                    n : in integer32;
                    ipvt : out Standard_Integer_Vectors.Vector;
                    rcond : out octo_double );

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

  procedure lufco ( a : in out OctoDobl_Complex_VecVecs.VecVec;
                    n : in integer32;
                    ipvt : out Standard_Integer_Vectors.Vector;
                    rcond : out octo_double );

  -- DESCRIPTION :
  --   lufco factors a complex matrix by gaussian elimination
  --   and estimates the condition of the matrix.
  --   The matrix is given as a vector of n columns.
  --   Each column contains the range 1..n.

  procedure lusolve ( a : in OctoDobl_Complex_Matrices.Matrix;
                      n : in integer32;
                      ipvt : in Standard_Integer_Vectors.Vector;
                      b : in out OctoDobl_Complex_Vectors.Vector );

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

  procedure lusolve ( a : in OctoDobl_Complex_VecVecs.VecVec;
                      n : in integer32;
                      ipvt : in Standard_Integer_Vectors.Vector;
                      b : in out OctoDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Version of lusolve for matrices defined as vectors of columns.

  procedure Triangulate ( a : in out OctoDobl_Complex_Matrices.Matrix;
                          tol : in double_float;
                          n,m : in integer32 );

  -- DESCRIPTION :
  --   Triangulate makes the n*m complex matrix a upper triangular using
  --   Gaussian elimination.

  -- ON ENTRY :
  --   a       a complex matrix(1..n,1..m);
  --   tol     tolerance to decide whether number is zero;
  --   n       the number of rows of a;
  --   m       the number of columns of a.

  -- ON RETURN :
  --   a       the triangulated matrix.

  procedure Diagonalize ( a : in out OctoDobl_Complex_Matrices.Matrix;
                          n,m : in integer32 );

  -- DESCRIPTION :
  --   Diagonalize makes the n*m complex matrix a diagonal using
  --   the Gauss-Jordan method.

  -- ON ENTRY :
  --   a       a complex matrix(1..n, 1..m);
  --   n       the number of rows of a;
  --   m       the number of columns of a.

  -- ON RETURN :
  --   a       the diagonalized matrix

-- TO TEST THE LU FACTORIZATION :

  function Permutation_Matrix
              ( ipvt : Standard_Integer_Vectors.Vector )
              return Standard_Natural_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the permutation matrix defined by the pivot selection in ipvt.

  function Permute ( P : Standard_Natural_Matrices.Matrix;
                     A : OctoDobl_Complex_Matrices.Matrix )
                   return OctoDobl_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the matrix P*A.

  function Lower_Diagonal ( A : OctoDobl_Complex_Matrices.Matrix )
                          return OctoDobl_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the lower diagonal part of A, with ones on the diagonal.

  function Upper_Diagonal ( A : OctoDobl_Complex_Matrices.Matrix )
                          return OctoDobl_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the upper diagonal part of A.

  procedure Permute_Lower
              ( L : in out OctoDobl_Complex_Matrices.Matrix;
                ipvt : in Standard_Integer_Vectors.Vector );

  -- DESCRIPTION :
  --   Applies the pivot selection to L = Lower_Diagonal(LU),
  --   the multipliers used in the LU factorization.
  --   This permutation is necessary to ensure that P*A = L*U.

end OctoDobl_Complex_Linear_Solvers;
