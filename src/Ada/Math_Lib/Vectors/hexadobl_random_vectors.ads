with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Hexa_Double_Vectors;
with HexaDobl_Complex_Vectors;

package HexaDobl_Random_Vectors is

-- DESCRIPTION :
--   Offers routines to generate vectors of random hexa double numbers.
--   Either the seed remains hidden for the user, or the user can manage
--   the seed for independent and/or reproducible sequences of numbers.

  function Random_Vector ( first,last : integer32 )
                         return Hexa_Double_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range first..last with random hexa doubles
  --   in the interval [-1,+1].

  procedure Random_Vector
              ( seed : in out integer32;
                v : out Hexa_Double_Vectors.Vector );

  -- DESRIPTION :
  --   Given a seed, generates a vector v of random hexa doubles,
  --   with each entry in [-1,+1].
  --   The seed is updated so the next call returns a different v.

  function Random_Vector ( first,last : integer32 )
                         return HexaDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range first..last with random hexa
  --   double complex numbers on the unit circle.

  procedure Random_Vector
              ( seed : in out integer32;
                v : out HexaDobl_Complex_Vectors.Vector );

  -- DESRIPTION :
  --   Given a seed, generates a vector v of random complex numbers of
  --   hexa double precision on the unit circle.
  --   The seed is updated so the next call returns a different v.

end HexaDobl_Random_Vectors;
