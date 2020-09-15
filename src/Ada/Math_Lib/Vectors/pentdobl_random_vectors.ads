with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Penta_Double_Vectors;
with PentDobl_Complex_Vectors;

package PentDobl_Random_Vectors is

-- DESCRIPTION :
--   Offers routines to generate vectors of random penta double numbers.
--   Either the seed remains hidden for the user, or the user can manage
--   the seed for independent and/or reproducible sequences of numbers.

  function Random_Vector ( first,last : integer32 )
                         return Penta_Double_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range first..last with random penta doubles
  --   in the interval [-1,+1].

  procedure Random_Vector
              ( seed : in out integer32;
                v : out Penta_Double_Vectors.Vector );

  -- DESRIPTION :
  --   Given a seed, generates a vector v of random penta doubles,
  --   with each entry in [-1,+1].
  --   The seed is updated so the next call returns a different v.

  function Random_Vector ( first,last : integer32 )
                         return PentDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range first..last with random penta
  --   double complex numbers on the unit circle.

  procedure Random_Vector
              ( seed : in out integer32;
                v : out PentDobl_Complex_Vectors.Vector );

  -- DESRIPTION :
  --   Given a seed, generates a vector v of random complex numbers of
  --   penta double precision on the unit circle.
  --   The seed is updated so the next call returns a different v.

end PentDobl_Random_Vectors;
