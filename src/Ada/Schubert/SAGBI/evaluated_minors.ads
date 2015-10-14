with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Floating_Matrices;
with Standard_Complex_Matrices;
with Brackets;                           use Brackets;

package Evaluated_Minors is

-- DESCRIPTION :
--   This package provides facilities to evaluate minors.

  function Determinant ( m : Standard_Floating_Matrices.Matrix )
                       return double_float;

  function Determinant ( m : Standard_Floating_Matrices.Matrix; b : Bracket )
                       return double_float;

  function Determinant ( m : Standard_Complex_Matrices.Matrix )
                       return Complex_Number;

  function Determinant ( m : Standard_Complex_Matrices.Matrix; b : Bracket )
                       return Complex_Number;

  -- DESCRIPTION :
  --   Computes the determinant of the sub-matrix of m, obtained by selecting
  --   those rows from the matrix m that are entries in the bracket b.
  --   If b is omitted, then the determinant of m is returned.

  -- REQUIRED :
  --   The matrix m should at least have as many columns as b'range and
  --   as many rows as b(b'last).  If b is omitted, then m must be square.

  function Determinant ( m : Standard_Complex_Matrices.Matrix;
                         rows,cols : Bracket ) return Complex_Number;

  -- DESCRIPTION :
  --   Returns the determinant of the matrix obtained from selecting
  --   those rows and columns as defined by the brackets rows and cols.

  -- REQUIRED : m'range(1) contains all indices of rows
  --   and m'range(2) contains all indices of cols.

end Evaluated_Minors;
