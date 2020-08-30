with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Triple_Double_Matrices;
with TripDobl_Complex_Matrices;

package TripDobl_Random_Matrices is

-- DESCRIPTION :
--   Offers routines to generate matrices of random triple doubles.

  function Random_Matrix ( n,m : natural32 )
                         return Triple_Double_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns a matrix of range 1..n,1..m with random triple doubles.

  function Random_Matrix ( n,m : natural32 )
                         return TripDobl_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns a matrix of range 1..n,1..m
  --   with random complex triple double numbers.

end TripDobl_Random_Matrices;
