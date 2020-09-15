with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Triple_Double_Vectors;
with TripDobl_Complex_Vectors;

package TripDobl_Random_Vectors is

-- DESCRIPTION :
--   Offers routines to generate vectors of random triple double numbers.
--   Either the seed remains hidden for the user, or the user can manage
--   the seed for independent and/or reproducible sequences of numbers.

  function Random_Vector ( first,last : integer32 )
                         return Triple_Double_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range first..last with random triple doubles
  --   in the interval [-1,+1].

  procedure Random_Vector
              ( seed : in out integer32;
                v : out Triple_Double_Vectors.Vector );

  -- DESRIPTION :
  --   Given a seed, generates a vector v of random triple doubles,
  --   with each entry in [-1,+1].
  --   The seed is updated so the next call returns a different v.

  function Random_Vector ( first,last : integer32 )
                         return TripDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range first..last with random triple
  --   double complex numbers on the unit circle.

  procedure Random_Vector
              ( seed : in out integer32;
                v : out TripDobl_Complex_Vectors.Vector );

  -- DESRIPTION :
  --   Given a seed, generates a vector v of random complex numbers of
  --   triple double precision on the unit circle.
  --   The seed is updated so the next call returns a different v.

end TripDobl_Random_Vectors;
