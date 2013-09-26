with Standard_Integer_Numbers;         use Standard_Integer_Numbers;
with Standard_Floating_Numbers;        use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_Matrices;

package Basis_Exchanges is

-- DESCRIPTION :
--   This package provides facilities for updating basis inverses
--   when rows of the basis are exchanged.
--   The update procedure is isolated in this package to enable massive
--   black-box testing when enumerating all possible bases.

  procedure Initial_Basis
                ( n,m : in integer32;
                  mat : in Standard_Floating_Matrices.Matrix;
                  tol : in double_float;
                  binv : out Standard_Floating_Matrices.Matrix;
                  cols : out Standard_Integer_Vectors.Vector;
                  fail : out boolean );

  -- DESCRIPTION :
  --   Generates an initial inverse of the basis in the matrix binv,
  --   that is the inverse to the selected columns of the given matrix.

  -- REQUIRED : binv'range(1) = binv'range(2) = 1..n;
  --   mat'range(1) contains 1..n, cols'range and mat'range(2) contains 1..m.
  
  -- ON ENTRY :
  --   n          number of variables;
  --   m          number of columns in the matrix;
  --   mat        column-oriented coefficient matrix;
  --   tol        to decide on zero of a real number.

  -- ON RETURN :
  --   binv       inverse of the basis to selected columns of the matrix mat;
  --   cols       selected columns from the given matrix mat;
  --   fail       when true then the matrix mat is not of full rank.

  function Solve ( n : integer32;
                   binv : Standard_Floating_Matrices.Matrix;
                   rhs : Standard_Floating_Vectors.Vector;
                   cols : Standard_Integer_Vectors.Vector )
                 return Standard_Floating_Vectors.Vector;

  -- DESCRIPTION :
  --   Uses the basis inverse binv to solve the system mat*x = rhs,
  --   for the given selection of the columns in mat and rhs.

  -- REQUIRED : rhs'range = mat'range(2).

  procedure Update ( n,m : in integer32;
                     binv : in out Standard_Floating_Matrices.Matrix;
                     mat : in Standard_Floating_Matrices.Matrix;
                     cols : in Standard_Integer_Vectors.Vector;
                     binv_row,mat_col : in integer32; tol : in double_float;
                     fail : out boolean );

  -- DESCRIPTION :
  --   Updates the inverse of the basis interchanging a row with a
  --   column from the coefficient matrix.
  --   The rows in binv are perpendicular to columns in the matrix mat.

  -- REQUIRED : binv'range(2) = binv'range(1) = 1..n;
  --   mat'range(1) contains 1..n, mat'range(2) contains 1..m.

  -- ON ENTRY :
  --   n          number of variables <= number of rows in mat;
  --   m          number of constraints <= number of columns in mat;
  --   binv       inverse of a basis taken from columns in mat;
  --   mat        column-oriented coefficient matrix;
  --   cols       current selection of columns in the coefficient matrix;
  --   binv_row   index to row in binv that has to leave the basis;
  --   mat_col    index to row in binv that comes into the basis;
  --   tol        everything smaller than tol is regarded as zero.

  -- ON RETURN :
  --   binv       inverse of basis with new row binv_row replaced
  --              by the row mat_col of the coefficient matrix;
  --   fail       no matching entry found in the row mat_col.

-- RESTARTING PROCEDURES :

  procedure Column_Basis
                ( n : in integer32;
                  mat : in Standard_Floating_Matrices.Matrix;
                  cols : in Standard_Integer_Vectors.Vector;
                  binv : out Standard_Floating_Matrices.Matrix;
                  fail : out boolean );

  -- DESCRIPTION :
  --   Computes the basis from a selection of columns of the matrix.
  --   This is useful to deal with the accumulation of rounding errors
  --   in the successive updates of the bases.

  -- REQUIRED : binv'range(1) = binv'range(2) = cols'range = 1..n;
  --   mat'range(1) contains 1..n and mat'range(2) contains entries of cols.
  
  -- ON ENTRY :
  --   n          number of variables;
  --   mat        column-oriented coefficient matrix;
  --   cols       selected columns from the given matrix mat.

  -- ON RETURN :
  --   binv       inverse of the basis to selected columns of the matrix mat;
  --   fail       if true then the selected columns form singular matrix.

  procedure Column_Solve
                ( n : in integer32;
                  mat : in Standard_Floating_Matrices.Matrix;
                  cols : in Standard_Integer_Vectors.Vector;
                  rhs : in Standard_Floating_Vectors.Vector;
                  sol : out Standard_Floating_Vectors.Vector;
                  fail : out boolean );

  -- DESCRIPTION :
  --   Computes the solution for the system mat*x = rhs with a given
  --   selection of columns in mat.  This is useful when the inverse of the
  --   basis is contaminated by too many accumulated rounding errors.

  -- REQUIRED : sol'range, cols'range and mat'range(1) must contain 1..n.

  -- ON ENTRY :
  --   n          number of variables;
  --   mat        column-oriented coefficient matrix;
  --   cols       selected columns from the given matrix mat;
  --   rhs        right-hand side vector of range mat'range(2).

  -- ON RETURN :
  --   sol        if not fail, then this is the solution to the system;
  --   fail       if true then the selected columns form a singular matrix.

-- ENUMERATOR :

  generic

    with procedure Report ( invbas : in Standard_Floating_Matrices.Matrix;
                            matcls : in Standard_Integer_Vectors.Vector;
                            level : in integer32; continue : out boolean );

    -- DESCRIPTION :
    --   This procedure is called after each update of the basis.

    -- ON ENTRY :
    --   invbas   newly updated inverse of the basis;
    --   matcls   columns of coefficient matrix used;
    --   level    deepness level of recursion in enumeration.

    -- ON RETURN :
    --   continue if true, the iteration continues, otherwise it stops.

  procedure Enumerate_Basis_Inverses
                ( n,m : in integer32;
                  mat,binv : in Standard_Floating_Matrices.Matrix;
                  mat_cols : in Standard_Integer_Vectors.Vector;
                  tol : in double_float );

  -- DESCRIPTION :
  --   Enumerates all possible inverses of bases taken from rows
  --   of the coefficient matrix.  After each update Report is called.

  -- REQUIRED : binv'range(2) = binv'range(1) = 1..n;
  --   mat'range(1) contains 1..n,
  --   mat'range(2) and mat_cols'range contain 1..m.

  -- ON ENTRY :
  --   mat        column-oriented coefficient matrix;
  --   binv       inverse of a basis from rows in the coefficient matrix;
  --   mat_cols   columns from the coefficient matrix that make up the basis;
  --   tol        tolerance for the absolute value of nonzero numbers.

end Basis_Exchanges;
