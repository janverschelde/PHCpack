with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Octo_Double_Matrices;
with OctoDobl_Complex_Matrices;

package OctoDobl_Random_Matrices is

-- DESCRIPTION :
--   Offers routines to generate matrices of random octo doubles.

  function Random_Matrix ( n,m : natural32 )
                         return Octo_Double_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns a matrix of range 1..n,1..m with random octo doubles.

  function Random_Matrix ( n,m : natural32 )
                         return OctoDobl_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns a matrix of range 1..n,1..m
  --   with random complex octo double numbers.

end OctoDobl_Random_Matrices;
