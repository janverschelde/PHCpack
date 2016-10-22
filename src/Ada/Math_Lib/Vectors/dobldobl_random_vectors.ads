with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Double_Double_Vectors;
with DoblDobl_Complex_Vectors;

package DoblDobl_Random_Vectors is

-- DESCRIPTION :
--   Offers routines to generate vectors of random double double numbers.
--   Either the seed remains hidden for the user, or the user can manage
--   the seed for independent and/or reproducible sequences of numbers.

  function Random_Vector ( first,last : integer32 )
                         return Double_Double_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range first..last with random double doubles
  --   in the interval [-1,+1].

  procedure Random_Vector
              ( seed : in out integer32;
                v : out Double_Double_Vectors.Vector );

  -- DESRIPTION :
  --   Given a seed, generates a vector v of random double doubles,
  --   with each entry in [-1,+1].
  --   The seed is updated so the next call returns a different v.

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

  procedure Random_Vector
              ( seed : in out integer32;
                v : out DoblDobl_Complex_Vectors.Vector );

  -- DESRIPTION :
  --   Given a seed, generates a vector v of random complex numbers
  --   of double double precision with modulus one.
  --   The seed is updated so the next call returns a different v.

  function Random_Vector ( first,last : integer32; m : natural32 )
                         return DoblDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range first..last with random 
  --   double double complex numbers with modulus in [10^(-m),10^(+m)].

end DoblDobl_Random_Vectors;
