with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Matrices;
with Standard_Integer64_Matrices;

package Standard_Smith_Normal_Form is

-- DESCRIPTION :
--   This package provides a diagonalization routine for a matrix of
--   standard integer numbers, along with some utilities.

  function Identity ( n : natural32 ) return Standard_Integer_Matrices.Matrix;
  function Identity ( n : natural32 ) return Standard_Integer64_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the n-by-n identity matrix.

  function Diagonal ( A : Standard_Integer_Matrices.Matrix ) return boolean;
  function Diagonal ( A : Standard_Integer64_Matrices.Matrix ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the matrix A is diagonal, returns false otherwise.

  function Rank_of_Diagonal_Matrix
             ( D : Standard_Integer_Matrices.matrix ) return natural32;
  function Rank_of_Diagonal_Matrix
             ( D : Standard_Integer64_Matrices.matrix ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of nonzero consecutive elements on
  --   the diagonal matrix D.

  -- REQUIRED : d is a diagonal matrix and all its nonzero elements
  --   on the diagonal are consecutive.

  procedure Diagonalize
              ( U,A,V : in out Standard_Integer_Matrices.Matrix );
  procedure Diagonalize
              ( U,A,V : in out Standard_Integer64_Matrices.Matrix );

  -- DESCRIPTION :
  --   Returns the Smith normal form of the matrix, A := U*A*V.

  -- ON ENTRY :
  --   U        identity matrix of range a'range(1);
  --   A        matrix we want to diagonalize;
  --   V        identity matrix of range a'range(2).

  -- ON RETURN :
  --   U        unimodular matrix of left multiplications;
  --   A        diagonal matrix, equals U*A*V;
  --   V        unimodular matrix of right multiplications.

end Standard_Smith_Normal_Form;
