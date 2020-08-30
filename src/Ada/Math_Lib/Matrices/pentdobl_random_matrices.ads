with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Penta_Double_Matrices;
with PentDobl_Complex_Matrices;

package PentDobl_Random_Matrices is

-- DESCRIPTION :
--   Offers routines to generate matrices of random penta doubles.

  function Random_Matrix ( n,m : natural32 )
                         return Penta_Double_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns a matrix of range 1..n,1..m with random penta doubles.

  function Random_Matrix ( n,m : natural32 )
                         return PentDobl_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns a matrix of range 1..n,1..m
  --   with random complex penta double numbers.

end PentDobl_Random_Matrices;
