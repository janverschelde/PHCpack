with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Deca_Double_Matrices;
with DecaDobl_Complex_Matrices;

package DecaDobl_Random_Matrices is

-- DESCRIPTION :
--   Offers routines to generate matrices of random deca doubles.

  function Random_Matrix ( n,m : natural32 )
                         return Deca_Double_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns a matrix of range 1..n,1..m with random deca doubles.

  function Random_Matrix ( n,m : natural32 )
                         return DecaDobl_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns a matrix of range 1..n,1..m
  --   with random complex deca double numbers.

end DecaDobl_Random_Matrices;
