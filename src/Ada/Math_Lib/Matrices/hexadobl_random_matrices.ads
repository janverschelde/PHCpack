with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Hexa_Double_Matrices;
with HexaDobl_Complex_Matrices;

package HexaDobl_Random_Matrices is

-- DESCRIPTION :
--   Offers routines to generate matrices of random hexa doubles.

  function Random_Matrix ( n,m : natural32 )
                         return Hexa_Double_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns a matrix of range 1..n,1..m with random hexa doubles.

  function Random_Matrix ( n,m : natural32 )
                         return HexaDobl_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns a matrix of range 1..n,1..m
  --   with random complex hexa double numbers.

end HexaDobl_Random_Matrices;
