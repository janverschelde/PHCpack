with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Double_Double_Vectors;
with DoblDobl_Complex_Vectors;

package DoblDobl_Random_Vectors is

-- DESCRIPTION :
--   Offers routines to generate vectors of random double double numbers.

  function Random_Vector ( first,last : integer32 )
                         return Double_Double_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range first..last with random double doubles
  --   in the interval [-1,+1].

  function Random_Vector ( first,last : integer32; m : natural32 )
                         return Double_Double_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range first..last with random double doubles
  --   with absolute value in [10^(-m),10^(+m)].

  function Random_Vector ( first,last : integer32 )
                         return DoblDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range first..last with random 
  --   double double complex numbers on the unit circle.

  function Random_Vector ( first,last : integer32; m : natural32 )
                         return DoblDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range first..last with random 
  --   double double complex numbers with modulus in [10^(-m),10^(+m)].

end DoblDobl_Random_Vectors;
