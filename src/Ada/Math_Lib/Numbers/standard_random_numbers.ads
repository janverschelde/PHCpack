with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;

package Standard_Random_Numbers is

-- DESCRIPTION :
--   This package provides random number generators for machine numbers.

  procedure Set_Seed ( n : natural32 );

  -- DESCRIPTION :
  --   Sets the seed to the number n.

  function Get_Seed return integer32;

  -- DESCRIPTION :
  --   Returns the seed used to generate random numbers.

  function Random ( lower,upper : integer32 ) return integer32;
  function Random ( lower,upper : integer64 ) return integer64;

  -- DESCRIPTION :
  --   lower <= Random(lower,upper) <= upper, randomly generated.

  function Random return double_float;

  -- DESCRIPTION :
  --   Returns a random floating-point number in [-1,1].

  function Random return Complex_Number;

  -- DESCRIPTION :
  --   Returns a random complex number, with real and imaginary
  --   part in [-1,1].

  function Random ( modulus : double_float ) return Complex_Number;

  -- DESCRIPTION :
  --   Generates a random complex number with a given modulus,
  --   so only the argument angle will be chosen at random.

  function Random1 return Complex_Number;

  -- DESCRIPTION :
  --   Generates a random complex number with modulus one.

  function Random_Magnitude ( m : natural32 ) return double_float;
 
  -- DESCRIPTION :
  --   Returns a random double float in [10^(-m),10^(+m)].

  function Random_Magnitude ( m : natural32 ) return Complex_Number;

  -- DESCRIPTION :
  --   Returns a random complex number with modulus in [10^(-m),10^(+m)].

end Standard_Random_Numbers;
