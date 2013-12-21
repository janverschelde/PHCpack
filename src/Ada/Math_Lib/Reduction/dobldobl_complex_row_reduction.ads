with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with DoblDobl_Complex_Matrices;          use DoblDobl_Complex_Matrices;

package DoblDobl_Complex_Row_Reduction is

-- DESCRIPTION :
--   This package offers routines to implement an incremental row reduction
--   on a matrix with standard complex numbers to determine whether the
--   matrix is singular or not.  This package was developed to solve
--   linear-product start systems efficiently.

  function Start_Pivots 
              ( n : integer32 ) return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the vector with entries 1 2 .. n to serve as start pivots
  --   in the beginning of the row reduction.

  function Pivot_in_Row ( A : Matrix; i,k : integer32;
                          tol : double_float ) return integer32;

  -- DESCRIPTION :
  --   Returns the pivot in row i of A, starting at the k-th element.

  -- ON ENTRY :
  --   A        a matrix of complex numbers, with defined i-th row;
  --   i        row to search for pivot;
  --   k        column index where to start the search;
  --   tol      if |x| < tol, then x is considered to be zero.
  -- 
  -- ON RETURN :
  --   Index of the first nonzero element on the i-th row of A,
  --   starting at column k; if the i-th row starting at column k
  --   is zero, then 0 is returned.

  procedure Swap_Columns
              ( A : in out Matrix; i,j : in integer32;
                piv : in out Standard_Integer_Vectors.Vector );

  -- DESCRIPTION :
  --   Swaps the elements on the i-th and j-th column on row i of A
  --   and accordingly in the pivots piv.

  procedure Divide_Row_by_Pivot ( A : in out Matrix; i : in integer32 );

  -- DESCRIPTION :
  --   Since A(i,i) /= 0, we divide all elements in the i-th row of A
  --   following the i-th column by A(i,i).

  procedure Eliminate ( A : in out Matrix; i : in integer32;
                        tol : in double_float );

  -- DESCRIPTION :
  --   The first i-1 rows of A are in upper triangular form 
  --   with ones on the diagonal.
  --   After Eliminate, the first i rows of A are triangular.

  procedure Reduce_Row 
              ( A : in out Matrix; i : in integer32;
                pivots : in out Standard_Integer_Vectors.Vector;
                tol : in double_float; singular : out boolean );
  procedure Reduce_Row 
              ( file : in file_type; A : in out Matrix; i : in integer32;
                pivots : in out Standard_Integer_Vectors.Vector;
                tol : in double_float; singular : out boolean );

  -- DESCRIPTION :
  --   Updates the i-th row of A to make it upper triangular.   

  -- REQUIRED :
  --   Rows 1 to i-1 of A are upper triangular with ones on the diagonal. 

  -- ON ENTRY :
  --   file     for intermediate output;
  --   A        matrix, already upper triangular in the first i-1 rows;
  --   i        current row to be updated;
  --   pivots   pivoting information from column interchanges;
  --   tol      if |x| < tol, then x is considered to be zero.

  -- ON RETURN :
  --   A        matrix upper triangular in the first i rows;
  --   singular is true if the i-th row after row reduction is zero,
  --            false otherwise.

end DoblDobl_Complex_Row_Reduction;
