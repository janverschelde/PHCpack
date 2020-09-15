with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Deca_Double_Vectors;
with DecaDobl_Complex_Vectors;

package DecaDobl_Random_Vectors is

-- DESCRIPTION :
--   Offers routines to generate vectors of random deca double numbers.
--   Either the seed remains hidden for the user, or the user can manage
--   the seed for independent and/or reproducible sequences of numbers.

  function Random_Vector ( first,last : integer32 )
                         return Deca_Double_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range first..last with random deca doubles
  --   in the interval [-1,+1].

  procedure Random_Vector
              ( seed : in out integer32;
                v : out Deca_Double_Vectors.Vector );

  -- DESRIPTION :
  --   Given a seed, generates a vector v of random deca doubles,
  --   with each entry in [-1,+1].
  --   The seed is updated so the next call returns a different v.

  function Random_Vector ( first,last : integer32 )
                         return DecaDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range first..last with random deca
  --   double complex numbers on the unit circle.

  procedure Random_Vector
              ( seed : in out integer32;
                v : out DecaDobl_Complex_Vectors.Vector );

  -- DESRIPTION :
  --   Given a seed, generates a vector v of random complex numbers of
  --   deca double precision on the unit circle.
  --   The seed is updated so the next call returns a different v.

end DecaDobl_Random_Vectors;
