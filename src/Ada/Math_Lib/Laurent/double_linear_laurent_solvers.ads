with text_io;                           use text_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_Matrices;
with Standard_Complex_VecVecs;
with Standard_Complex_VecVecVecs;

package Double_Linear_Laurent_Solvers is

-- DESCRIPTION :
--   Linear systems defined by coefficient matrices and right hand side
--   vectors of Laurent series can be solved by an LU factorization.
--   All computations happen in double precision.

  procedure Allocate_Series_Coefficients
              ( dim,deg : in integer32;
                cff : out Standard_Complex_VecVecs.Link_to_VecVec );

  -- DESCRIPTION :
  --   Returns in cff the coefficients of dim series of degree deg,
  --   all equal to zero.

  procedure Write ( e : in Standard_Integer_Matrices.Matrix;
                    c : in Standard_Complex_VecVecVecs.Link_to_VecVecVec;
                    s : in string := "A" );

  -- DESCRIPTION :
  --   Writes the matrix of Laurent series, defined by
  --   the leading exponents in e and coefficients in c.
  --   The string is used as the name of the matrix.

  -- REQUIRED :
  --   e'range(1) = c'range(1) and e'range(2) = c'range(2).

  procedure Write ( e : in Standard_Integer_Vectors.Vector;
                    c : in Standard_Complex_VecVecs.Link_to_VecVec;
                    s : in string := "v" );
  procedure Write ( file : in file_type;
                    e : in Standard_Integer_Vectors.Vector;
                    c : in Standard_Complex_VecVecs.Link_to_VecVec;
                    s : in string := "v" );

  -- DESCRIPTION :
  --   Writes the vector of Laurent series to standard output,
  --   or to file.  The vector of Laurent series is defined by
  --   the leading exponents in e and coefficients in c.
  --   The string s is used as the name of the vector.

  -- REQUIRED : e'range = c'range.

  procedure Matrix_Vector_Product
              ( d : in integer32;
                eA : in Standard_Integer_Matrices.Matrix;
                cA : in Standard_Complex_VecVecVecs.Link_to_VecVecVec;
                ex : in Standard_Integer_Vectors.Vector;
                cx : in Standard_Complex_VecVecs.Link_to_VecVec;
                ey : out Standard_Integer_Vectors.Vector;
                cy : out Standard_Complex_VecVecs.Link_to_VecVec );

  -- DESCRIPTION :
  --   Returns the product of a matrix of Laurent series
  --   with a vector of Laurent series.

  -- REQUIRED :
  --   eA'range(1) = ey'range and eA'range(2) = ex'range,
  --   cA'range(1) = cy'range and cA'range(2) = cx'range.

  -- ON ENTRY :
  --   d        only coefficients in the range 0 to d are considered;
  --   eA       leading exponents of the series in the matrix A;
  --   cA       coefficients of the series in the matrix A;
  --   ex       leading exponents of the series in the vector x;
  --   cx       coefficients of the series in the vector x;
  --   cy       space allocated for the coefficients of the product.

  -- ON RETURN :
  --   ey       leading exponents of the series of the product;
  --   cy       leading coefficients of the series of the product.

  procedure Forward_Substitution
              ( d : in integer32;
                eL : in Standard_Integer_Matrices.Matrix;
                cL : in Standard_Complex_VecVecVecs.Link_to_VecVecVec;
                eb : in Standard_Integer_Vectors.Vector;
                cb : in Standard_Complex_VecVecs.Link_to_VecVec;
                ex : out Standard_Integer_Vectors.Vector;
                cx : out Standard_Complex_VecVecs.Link_to_VecVec );

  -- DESCRIPTION :
  --   Applies forward substitution to solve a lower triangular system
  --   with ones on the diagonal.

  -- REQUIRED :
  --   The matrix is square and all ranges are compatible.

  -- ON ENTRY :
  --   d        only coefficients in the range 0 to d are considered;
  --   eL       leading exponents in the lower triangular matrix L;
  --   cL       coefficients of the series in the matrix L;
  --   eb       leading exponents of the right hand side vector b;
  --   cb       coefficients of the series in the vector b;
  --   cx       space allocated for the coefficients of the solution.

  -- ON RETURN :
  --   ex       leading exponents of the series of the solution;
  --   cx       leading coefficients of the series of the solution.

  procedure Backward_Substitution
              ( d : in integer32;
                eU : in Standard_Integer_Matrices.Matrix;
                cU : in Standard_Complex_VecVecVecs.Link_to_VecVecVec;
                eb : in Standard_Integer_Vectors.Vector;
                cb : in Standard_Complex_VecVecs.Link_to_VecVec;
                ex : out Standard_Integer_Vectors.Vector;
                cx : out Standard_Complex_VecVecs.Link_to_VecVec );

  -- DESCRIPTION :
  --   Applies forward substitution to solve an upper triangular system
  --   with general, nonzero elements on the diagonal.

  -- REQUIRED :
  --   The matrix is square and all ranges are compatible.

  -- ON ENTRY :
  --   d        only coefficients in the range 0 to d are considered;
  --   eU       leading exponents in the upper triangular matrix U;
  --   cU       coefficients of the series in the matrix U;
  --   eb       leading exponents of the right hand side vector b;
  --   cb       coefficients of the series in the vector b;
  --   cx       space allocated for the coefficients of the solution.

  -- ON RETURN :
  --   ex       leading exponents of the series of the solution;
  --   cx       leading coefficients of the series of the solution.

  function Pivot_Row
              ( nrows,row : in integer32;
                lead : in Standard_Integer_Matrices.Matrix;
                column : in Standard_Complex_VecVecs.Link_to_VecVec )
              return integer32;

  -- DESCRIPTION :
  --   Returns the row of the pivot in the current column.
  --   The pivot is determined first by the smallest leading exponent.
  --   If the leading exponents agree, then the largest coefficient
  --   determines the pivot row.

  -- ON ENTRY :
  --   nrows    number of rows;
  --   row      index of the current row;
  --   lead     leading exponents of the series;
  --   column   coefficients of the series in the current column.

  procedure Swap_Rows
              ( ncols,row,pivrow : in integer32;
                lead : in out Standard_Integer_Matrices.Matrix;
                cffs : in Standard_Complex_VecVecVecs.Link_to_VecVecVec;
                pivots : in out Standard_Integer_Vectors.Vector );

  -- DESCRIPTION :
  --   Swaps the rows, as defined by the index of the pivot row.

  -- ON ENTRY :
  --   ncols    number of columns;
  --   row      index of the current row;
  --   pivrow   index of the pivot row, different from row;
  --   lead     leading exponents of the series;
  --   cffs     coefficients of the series;
  --   pivots   current values of the pivots.

  procedure LU_Factorization
              ( nrows,ncols,deg : in integer32;
                Alead : in out Standard_Integer_Matrices.Matrix;
                Acffs : in Standard_Complex_VecVecVecs.Link_to_VecVecVec;
                pivots : out Standard_Integer_Vectors.Vector );

  -- DESCRIPTION :
  --   An inplace LU factorization with pivoting.

  -- ON ENTRY :
  --   nrows    number of rows of the matrix;
  --   ncols    number of columns of the matrix;
  --   deg      degree of the series in the matrix;
  --   Alead    leading exponents of the series;
  --   Acffs    coefficients of the series;
  --   pivots   space for the pivots.

  -- ON RETURN :
  --   Alead    contains the leading exponents of the factors;
  --   Acffs    the coefficients of the factors L and U;
  --   pivots   are the pivots used.

end Double_Linear_Laurent_Solvers;
