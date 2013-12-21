with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Natural_Matrices;
with Abstract_Ring,Abstract_Ring.Field;
with Generic_Vectors,Generic_Matrices;

generic

  with package Ring is new Abstract_Ring(<>);
  with package Field is new Ring.Field(<>);
  with package Vectors is new Generic_Vectors(Ring);
  with package Matrices is new Generic_Matrices(Ring,Vectors);

package Generic_Floating_Linear_Solvers is

-- DESCRIPTION :
--   This package offers a few routines to solve linear systems of equations.
--   The code for lufac, lufco and lusolve is a literal translation from the
--   f77-linpack code.  We distinguish between static triangulators, which
--   take the matrix as a whole at once, and dynamic triangulators, which
--   allow to triangulate row after row.
 
  use Ring,Matrices;

-- STATIC TRIANGULATORS :

  procedure lufac ( a : in out Matrix; n : in integer32;
                    ipvt : out Standard_Integer_Vectors.Vector;
                    info : out integer32 );

  -- DESCRIPTION :
  --   lufac factors a float matrix by gaussian elimination

  --   lufac is usually called by lufco, but it can be called
  --   directly with a saving of time if rcond is not needed.
  --   (time for lufco) = (1 + 9/n)*(time for lufac).

  -- ON ENTRY :
  --   a            a floating point matrix(1..n,1..n) to be factored
  --   n            the dimension of the matrix a

  -- ON RETURN :
  --   a            an upper triangular matrix and the multipliers
  --                which were used to obtain it.
  --                The factorization can be written a = l*u where
  --                l is a product of permutation and unit lower
  --                triangular matrices and u is upper triangular.
  --   ipvt         a natural vector of pivot indices
  --   info         = 0  normal value
  --                = k  if u(k,k) = 0.0.
  --                     This is not an error for this routine,
  --                     but it does indicate that lusolve will
  --                     divide by zero if called.  Use rcond in
  --                     lufco for a reliable indication of singularity.

  procedure lufco ( a : in out Matrix; n : in integer32;
                    ipvt : out Standard_Integer_Vectors.Vector;
                    rcond : out number );

  -- DESCRIPTION :
  --   lufco factors a floating point matrix by gaussian elimination
  --   and estimates the condition of the matrix.

  -- If rcond is not needed, lufac is slightly faster.
  -- To solve a*x = b, follow lufco by lusolve.

  -- ON ENTRY :
  --   a           a floating point matrix(1..n,1..n) to be factored;
  --   n           the dimension of the matrix a.

  -- ON RETURN :
  --   a           an upper triangular matrix and the multipliers 
  --               which are used to obtain it.
  --               The factorization can be written a = l*u, where
  --               l is a product of permutation and unit lower triangular
  --               matrices and u is upper triangular.
  --   ipvt        a natural vector of pivot indices
  --   rcond       an estimate of the reciprocal condition of a.
  --               For the system a*x = b, relative perturbations
  --               in a and b of size epsilon may cause relative
  --               perturbations in x of size epsilon/rcond.
  --               If rcond is so small that the logical expression
  --                      1.0 + rcond = 1.0 
  --               is true, than a may be singular to working precision.
  --               In particular, rcond is zero if exact singularity is
  --               detected or the estimate underflows.

  procedure lusolve ( a : in Matrix; n : in integer32;
                      ipvt : in Standard_Integer_Vectors.Vector;
                      b : in out Vectors.Vector );

  -- DESCRIPTION :
  --   lusolve solves the system a*x = b using the factors
  --   computed by lufac or lufco

  -- ON ENTRY :
  --   a           a floating-point matrix(1..n,1..n), the output from
  --               lufac or lufco;
  --   n           the dimension of the matrix a;
  --   ipvt        the pivot vector from lufac or lufco;
  --   b           the right hand side vector.

  -- ON RETURN :
  --   b           the solution vector x.

  procedure Triangulate ( a : in out Matrix; n,m : in integer32 );

  -- DESCRIPTION :
  --   triangulate makes the n*m matrix a triangular using
  --   Gaussian elimination.

  -- ON ENTRY :
  --   a           a floating-point matrix(1..n,1..m);
  --   n           the number of rows of a;
  --   m           the number of columns of a.

  -- ON RETURN :
  --   a           the triangulated matrix.

  procedure Diagonalize ( a : in out Matrix; n,m : in integer32 );

  -- DESCRIPTION :
  --   Diagonalize makes the n*m floating-point matrix a diagonal.

  -- ON ENTRY :
  --   a           a floating-point matrix(1..n,1..m);
  --   n           the number of rows of a;
  --   m           the number of columns of a.

  -- ON RETURN :
  --   a           the diagonalized matrix.

-- DYNAMIC TRIANGULATORS :

  procedure Upper_Triangulate
                 ( row : in integer32; mat : in out Matrix; tol : in number;
                   ipvt : in out Standard_Integer_Vectors.Vector;
                   pivot : out integer32 );

  -- DESCRIPTION :
  --   Makes the matrix upper triangular by updating the current row of
  --   the matrix.  If pivot = 0 on return, then the matrix is singular.

  -- REQUIRED :
  --   The matrix is upper triangular up to current row, which means that
  --   abs(max(i,i)) > tol, for i in mat'first(1)..row-1.
  --   The matrix might have more columns than rows, but ipvt'range should
  --   match the range of the first columns in the matrix.

  procedure Upper_Triangulate
                 ( roweli : in integer32; elim : in Matrix; tol : in number;
                   rowmat : in integer32; mat : in out Matrix );
  procedure Upper_Triangulate
                 ( roweli : in integer32; elim : in Matrix; tol : in number;
                   firstrow,lastrow : in integer32; mat : in out Matrix );

  -- DESCRIPTION :
  --   Using the matrix elim, the unknown roweli is eliminated in mat.

  -- ON ENTRY :
  --   roweli      current unknown to be eliminated;
  --   elim        elim(1..roweli,m'range(2)) is upper triangular;
  --   firstrow    indicates start in range of mat to be updated;
  --   lastrow     indicates end in range of mat to be updated;
  --   roweli      indicates the current unknown to be updated;
  --   mat         the unknows before roweli are already eliminated.

  -- ON RETURN :
  --   mat       the updated matrix after elimination of roweli.

  procedure Switch ( ipvt : in Standard_Integer_Vectors.Vector;
                     row : in integer32; mat : in out Matrix );

  -- DESCRIPTION :
  --   Applies the pivoting information ipvt to the row of the matrix.

  procedure Switch ( k,pivot,first,last : in integer32; mat : in out Matrix );

  -- DESCRIPTION :
  --   Swaps the column k with the pivot column for row first..last in mat.

  function Solve ( mat : Matrix; tol : number;
                   ipvt : Standard_Integer_Vectors.Vector )
                 return Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a solution to the homogenous linear system with inequalities
  --   in the rows of the matrix.

  -- REQUIRED : the matrix is upper triangular.

  function Solve ( n,col : integer32; mat : Matrix )
                 return Vectors.Vector;

  -- DESCRIPTION :
  --   Solves the system with as right-hand side the column indicated
  --   by the parameter col.

  -- REQUIRED :
  --   The matrix must be upper triangular.

-- TO TEST THE LU FACTORIZATION :

  function Permutation_Matrix
              ( ipvt : Standard_Integer_Vectors.Vector )
              return Standard_Natural_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the permutation matrix defined by the pivot selection in ipvt.

  function Permute ( P : Standard_Natural_Matrices.Matrix;
                     A : Matrix ) return Matrix;

  -- DESCRIPTION :
  --   Returns the matrix P*A.

  function Lower_Diagonal ( A : Matrix ) return Matrix;

  -- DESCRIPTION :
  --   Returns the lower diagonal part of A, with ones on the diagonal.

  function Upper_Diagonal ( A : Matrix ) return Matrix;

  -- DESCRIPTION :
  --   Returns the upper diagonal part of A.

  procedure Permute_Lower
              ( L : in out Matrix;
                ipvt : in Standard_Integer_Vectors.Vector );

  -- DESCRIPTION :
  --   Applies the pivot selection to L = Lower_Diagonal(LU),
  --   the multipliers used in the LU factorization.
  --   This permutation is necessary to ensure that P*A = L*U.

end Generic_Floating_Linear_Solvers;
