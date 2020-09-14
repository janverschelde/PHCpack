with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Triple_Double_Numbers;              use Triple_Double_Numbers;
with TripDobl_Complex_Numbers;           use TripDobl_Complex_Numbers;

package TripDobl_Random_Numbers is

-- DESCRIPTION :
--   This package provides random number generators for 
--   triple doubles and complex triple doubles.
--   The seed is set in the package Standard_Random_Numbers.
--   Alternatively, the user can manage the seed for independent
--   and/or reproducible sequences of random numbers.

  function Random return triple_double;

  -- DESCRIPTION :
  --   Returns a random triple double in [-1,+1].

  procedure Random_Triple_Double
              ( seed : in out integer32; f : out triple_double );

  -- DESCRIPTION :
  --   Given a seed, returns a random triple double f in [-1,+1].
  --   The seed is updated so the next call returns a different f.
 
  function Random return Complex_Number;

  -- DESCRIPTION :
  --   Returns a random complex triple double with
  --   real and imaginary parts in [-1,+1].

  procedure Random_Complex_Number
              ( seed : in out integer32; c : out Complex_Number );

  -- DESCRIPTION :
  --   Given a seed, returns a random complex number c, with real
  --   and imaginary parts both in [-1,+1].
  --   The seed is updated so the next call returns a different c.

  function Random1 return Complex_Number;

  -- DESCRIPTION :
  --   Returns a random complex triple double with modulus equal to one.

  procedure Random1_Complex_Number
              ( seed : in out integer32; c : out Complex_Number );

  -- DESCRIPTION :
  --   Given a seed, returns a random complex number c,
  --   with modulus equal to one.
  --   The seed is updated so the next call returns a different c.

end TripDobl_Random_Numbers; 
