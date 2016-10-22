with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with QuadDobl_Complex_Numbers;           use QuadDobl_Complex_Numbers;

package QuadDobl_Random_Numbers is

-- DESCRIPTION :
--   This package provides random number generators for 
--   quad doubles and complex quad doubles.
--   The seed is set in the package Standard_Random_Numbers.
--   Alternatively, the user can manage the seed for independent
--   and/or reproducible sequences of random numbers.

  function Random return quad_double;

  -- DESCRIPTION :
  --   Returns a random quad double in [-1,+1].

  procedure Random_Quad_Double
              ( seed : in out integer32; f : out quad_double );

  -- DESCRIPTION :
  --   Given a seed, returns a random quad double f in [-1,+1].
  --   The seed is updated so the next call returns a different f.
 
  function Random_Magnitude ( m : natural32 ) return quad_double;

  -- DESCRIPTION :
  --   Returns a random quad double with absolute value
  --   in [10^(-m),10^(+m)].

  function Random return Complex_Number;

  -- DESCRIPTION :
  --   Returns a random complex quad double with
  --   real and imaginary parts in [-1,+1].

  procedure Random_Complex_Number
              ( seed : in out integer32; c : out Complex_Number );

  -- DESCRIPTION :
  --   Given a seed, returns a random complex number c, with real
  --   and imaginary parts both in [-1,+1].
  --   The seed is updated so the next call returns a different c.

  function Random1 return Complex_Number;

  -- DESCRIPTION :
  --   Returns a random complex quad double
  --   with modulus equal to one.

  procedure Random1_Complex_Number
              ( seed : in out integer32; c : out Complex_Number );

  -- DESCRIPTION :
  --   Given a seed, returns a random complex number c,
  --   with modulus equal to one.
  --   The seed is updated so the next call returns a different c.

  function Random_Magnitude ( m : natural32 ) return Complex_Number;

  -- DESCRIPTION :
  --   Returns a random complex quad double of modulus in
  --   the interval [10^(-m),10^(+m)].

end QuadDobl_Random_Numbers;
