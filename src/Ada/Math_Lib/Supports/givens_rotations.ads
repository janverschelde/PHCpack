with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;          use Standard_Floating_Vectors;
with Standard_Floating_Matrices;         use Standard_Floating_Matrices;

package Givens_Rotations is

-- DESCRIPTION :
--   This package contains operations to construct and perform
--   Given rotations on vectors and matrices.
--   The procedures in this package are listed in order like
--   they have to be applied to write a right hand side vector
--   as a linear combination of the column vectors of a matrix.

  procedure Givens_Factors ( v : in Vector; i,j : in integer32;
                             cosi,sine : out double_float );
  -- DESCRIPTION :
  --   Computes the cosine and sine of the angle to be used in the
  --   Givens rotation on the vector.

  -- REQUIRED : v(i) /= 0 and v(j) /= 0, with i<j.

  procedure Givens_Factors ( mat : in Matrix; i,j : in integer32;
                             cosi,sine : out double_float );
  -- DESCRIPTION :
  --   Computes the cosine and sine of the angle to be used in the
  --   Givens rotation on the matrix.

  -- REQUIRED : mat(i,i) /= 0 and mat(j,i) /= 0, with i<j.

  procedure Givens_Rotation ( v : in out Vector; i,j : in integer32;
                              cosi,sine : in double_float );
  -- DESCRIPTION :
  --   Performs one Givens rotation on the vector
  --   with given cosine and sine of the angle.

  -- REQUIRED : v(i) /= 0 and v(j) /= 0, with i<j.

  procedure Givens_Rotation ( mat : in out Matrix; lastcol,i,j : in integer32;
                              cosi,sine : in double_float );
  procedure Givens_Rotation ( mat : in out Matrix; i,j : in integer32;
                              cosi,sine : in double_float );

  -- DESCRIPTION :
  --   Performs one Givens rotation on the matrix,
  --   with given cosine and sine of the angle.
  --   The parameter lastcol is the index of the last column of interest
  --   in the matrix.

  -- REQUIRED :
  --   mat(i,i) /= 0 and mat(j,i) /= 0, with i<j.
  --   It will be assumed that mat is already upper triangular up
  --   to the ith row and colum, i.e. mat(k,l) = 0, l<i, k>l.

  procedure Givens_Rotation ( v : in out Vector; i,j : in integer32 );

  -- DESCRIPTION :
  --   Performs one Givens rotation on the vector v.

  -- REQUIRED : v(i) /= 0 and v(j) /= 0, with i<j.

  procedure Givens_Rotation ( mat : in out Matrix; i,j : in integer32 );
  procedure Givens_Rotation ( mat : in out Matrix; lastcol,i,j : in integer32 );

  -- DESCRIPTION :
  --   Performs one Givens rotation on the matrix.
  --   The parameter lastcol is the index of the last column of interest
  --   in the matrix.

  -- REQUIRED :
  --   mat(i,i) /= 0 and mat(i,j) /= 0, with i<j.
  --   It will be assumed that mat is already upper triangular up
  --   to the ith row and colum, i.e. mat(k,l) = 0, l<i, k>l.

  procedure Upper_Triangulate
               ( row : in integer32; mat : in out Matrix;
                 tol : in double_float;
                 ipvt : in out Standard_Integer_Vectors.Vector;
                 pivot : out integer32 );

  -- DESCRIPTION :
  --   Makes the matrix upper triangular by updating the current row of
  --   the matrix.  If pivot = 0 on return, then the matrix is singular.

  -- REQUIRED :
  --   The matrix is upper triangular up to current row, which means that
  --   abs(max(i,i)) > tol, for i in mat'first(1)..row-1.

  procedure Upper_Triangulate ( mat : in out Matrix; tol : in double_float );
  procedure Upper_Triangulate ( mat : in out Matrix; tol : in double_float;
                                ipvt : out Standard_Integer_Vectors.Vector );
  -- DESCRIPTION :
  --   Makes the matrix upper triangular by means of Givens rotations.
  --   The parameter tol is used to determine whether an element is zero.
  --   Column interchangements are performed when necessary.
  --   The pivoting information will be returned in the vector ipvt.
  --   Stops when a zero row is encountered.

  -- REQUIRED : mat'first(1) = mat'first(2).

  procedure Upper_Triangulate ( mat : in out Matrix; rhs : in out Vector;
                                tol : in double_float );
  procedure Upper_Triangulate ( mat : in out Matrix; rhs : in out Vector;
                                tol : in double_float; 
                                ipvt : out Standard_Integer_Vectors.Vector );
  -- DESCRIPTION :
  --   Makes the matrix upper triangular and
  --   performs all Givens rotations as well on the right hand side vector.
  --   The parameter tol is used to determine whether an element is zero.
  --   Column interchangements are performed when necessary.
  --   The pivoting information will be returned in the vector ipvt.
  --   Stops when a zero row is encountered.

  -- REQUIRED : mat'first(1) = mat(first(2).

  procedure Solve ( mat : in Matrix; rhs : in Vector; tol : in double_float;
                    x : out Vector );

  -- DESCRIPTION :
  --   Solves the system defined by mat*x = rhs.
  --   In case rank(mat) < mat'length(1), then only 
  --   x(mat'first(1)..mat'first(1)+rank(mat)-1) will be meaningful.
  
  -- REQUIRED :
  --   mat is upper triangular and mat(i,i) /= 0, for i in 1..Rank(mat).

  -- NOTE :
  --   Eventually, when column interchangements where applied during
  --   the triangulation of the matrix, then x is not the decomposition
  --   of the right hand side vector in terms of the first column of
  --   the matrix.

  procedure Solve ( mat : in Matrix; rhs : in Vector; tol : in double_float;
                    rank : in integer32; x : out Vector );

  -- DESCRIPTION :
  --   Solves the system defined by mat*x = rhs.
  --   Here only the first rows in columns up to the given rank are
  --   considered.  So only x(mat'first(1)..mat'first(1)+rank-1) will
  --   be meaningful.

  -- REQUIRED :
  --   The matrix is upper triangular up to its rank
  --   and mat(i,i) /= 0, for i in 1..rank.

end Givens_Rotations;
