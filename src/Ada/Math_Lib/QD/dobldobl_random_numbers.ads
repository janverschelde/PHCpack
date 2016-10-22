with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with DoblDobl_Complex_Numbers;           use DoblDobl_Complex_Numbers;

package DoblDobl_Random_Numbers is

-- DESCRIPTION :
--   This package provides random number generators for 
--   double doubles and complex double doubles.
--   The seed is set in the package Standard_Random_Numbers.
--   Alternatively, the user can manage the seed for independent
--   and/or reproducible sequences of random numbers.

  function Random return double_double;

  -- DESCRIPTION :
  --   Returns a random double double in [-1,+1].

  procedure Random_Double_Double
              ( seed : in out integer32; f : out double_double );

  -- DESCRIPTION :
  --   Given a seed, returns a random double double f in [-1,+1].
  --   The seed is updated so the next call returns a different f.
 
  function Random_Magnitude ( m : natural32 ) return double_double;

  -- DESCRIPTION :
  --   Returns a random double double with absolute value
  --   in [10^(-m),10^(+m)].

  function Random return Complex_Number;

  -- DESCRIPTION :
  --   Returns a random complex double double with
  --   real and imaginary parts in [-1,+1].

  procedure Random_Complex_Number
              ( seed : in out integer32; c : out Complex_Number );

  -- DESCRIPTION :
  --   Given a seed, returns a random complex number c, with real
  --   and imaginary parts both in [-1,+1].
  --   The seed is updated so the next call returns a different c.

  function Random1 return Complex_Number;

  -- DESCRIPTION :
  --   Returns a random complex double double
  --   with modulus equal to one.

  procedure Random1_Complex_Number
              ( seed : in out integer32; c : out Complex_Number );

  -- DESCRIPTION :
  --   Given a seed, returns a random complex number c,
  --   with modulus equal to one.
  --   The seed is updated so the next call returns a different c.
 
  function Random_Magnitude ( m : natural32 ) return Complex_Number;

  -- DESCRIPTION :
  --   Returns a random complex double double of modulus in
  --   the interval [10^(-m),10^(+m)].

end DoblDobl_Random_Numbers;
