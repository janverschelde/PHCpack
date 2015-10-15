with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Floating_Matrices;
with Standard_Complex_Matrices;
with Double_Double_Matrices;
with DoblDobl_Complex_Matrices;
with Quad_Double_Matrices;
with QuadDobl_Complex_Matrices;
with Brackets;                           use Brackets;

package Evaluated_Minors is

-- DESCRIPTION :
--   This package provides facilities to evaluate minors.

  function Determinant ( m : Standard_Floating_Matrices.Matrix )
                       return double_float;
  function Determinant ( m : Double_Double_Matrices.Matrix )
                       return double_double;
  function Determinant ( m : Quad_Double_Matrices.Matrix )
                       return quad_double;

  -- DESCRIPTION :
  --   Returns the determinant of the real matrix m, 
  --   computed in double, double double, or quad double precision.

  -- REQUIRED : m'range(1) = m'range(2), the matrix is square.

  function Determinant ( m : Standard_Floating_Matrices.Matrix; b : Bracket )
                       return double_float;
  function Determinant ( m : Double_Double_Matrices.Matrix; b : Bracket )
                       return double_double;
  function Determinant ( m : Quad_Double_Matrices.Matrix; b : Bracket )
                       return quad_double;

  -- DESCRIPTION :
  --   Computes the determinant of the submatrix of m, obtained by selecting
  --   those rows from the matrix m that are entries in the bracket b,
  --   in double, double double, or quad double precision.

  -- REQUIRED :
  --   The matrix m should at least have as many columns as b'last and
  --   as many rows as b(b'last).

  function Determinant ( m : Standard_Complex_Matrices.Matrix )
                       return Standard_Complex_Numbers.Complex_Number;
  function Determinant ( m : DoblDobl_Complex_Matrices.Matrix )
                       return DoblDobl_Complex_Numbers.Complex_Number;
  function Determinant ( m : QuadDobl_Complex_Matrices.Matrix )
                       return QuadDobl_Complex_Numbers.Complex_Number;

  -- DESCRIPTION :
  --   Returns the determinant of the complex matrix m, 
  --   computed in double, double double, or quad double precision.

  -- REQUIRED : m'range(1) = m'range(2), the matrix is square.

  function Determinant ( m : Standard_Complex_Matrices.Matrix; b : Bracket )
                       return Standard_Complex_Numbers.Complex_Number;
  function Determinant ( m : DoblDobl_Complex_Matrices.Matrix; b : Bracket )
                       return DoblDobl_Complex_Numbers.Complex_Number;
  function Determinant ( m : QuadDobl_Complex_Matrices.Matrix; b : Bracket )
                       return QuadDobl_Complex_Numbers.Complex_Number;

  -- DESCRIPTION :
  --   Computes the determinant of the sub-matrix of m, obtained by selecting
  --   those rows from the matrix m that are entries in the bracket b,
  --   in double, double double, or quad double precision.

  -- REQUIRED :
  --   The matrix m should at least have as many columns as b'last and
  --   as many rows as b(b'last).

  function Determinant ( m : Standard_Complex_Matrices.Matrix;
                         rows,cols : Bracket )
                       return Standard_Complex_Numbers.Complex_Number;
  function Determinant ( m : DoblDobl_Complex_Matrices.Matrix;
                         rows,cols : Bracket )
                       return DoblDobl_Complex_Numbers.Complex_Number;
  function Determinant ( m : QuadDobl_Complex_Matrices.Matrix;
                         rows,cols : Bracket )
                       return QuadDobl_Complex_Numbers.Complex_Number;

  -- DESCRIPTION :
  --   Returns the determinant of the matrix obtained from selecting
  --   those rows and columns as defined by the brackets rows and cols,
  --   in double, double double, and quad double precision.

  -- REQUIRED : m'range(1) contains all indices of rows
  --   and m'range(2) contains all indices of cols.

end Evaluated_Minors;
