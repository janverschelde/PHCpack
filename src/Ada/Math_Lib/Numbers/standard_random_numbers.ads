with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;

package Standard_Random_Numbers is

-- DESCRIPTION :
--   This package provides random number generators for machine numbers.
--   Unless the seed is specified explicitly, the seed of the random
--   number generators are based on the process identification number.

  procedure Set_Seed ( n : natural32 );

  -- DESCRIPTION :
  --   Sets the seed to the number n.

  function Get_Seed return integer32;

  -- DESCRIPTION :
  --   Returns the seed used to generate random numbers.

  function Random ( probability : double_float := 0.5 ) return boolean;

  -- DESCRIPTION :
  --   Returns true with given probability.
  --   By default, the likelihood of true is equal to that of false. 
  --   With higher probability, true is more likely,
  --   with lower probability, false is more likely.

  function Random ( lower,upper : integer32 ) return integer32;
  function Random ( lower,upper : integer64 ) return integer64;

  -- DESCRIPTION :
  --   lower <= Random(lower,upper) <= upper, randomly generated.

  procedure Random_Integer_Number
              ( seed : in out integer32; lower,upper : in integer32;
                i : out integer32 );
  procedure Random_Integer_Number
              ( seed : in out integer32; lower,upper : in integer64;
                i : out integer64 );

  -- DESCRIPTION :
  --   Given a seed, returns a random integer i in [lower, upper].
  --   The seed is updated so the next call returns another i.

  function Random return double_float;

  -- DESCRIPTION :
  --   Returns a random floating-point number in [-1,1].

  procedure Random_Double_Float 
              ( seed : in out integer32; f : out double_float );

  -- DESCRIPTION :
  --   Returns a random floating-point number f in [-1,1],
  --   using the given seed.  The seed is updated so the next
  --   call to this function will return another number in f.

  function Random return Complex_Number;

  -- DESCRIPTION :
  --   Returns a random complex number, with real and imaginary
  --   parts both in [-1,1].

  procedure Random_Complex_Number
              ( seed : in out integer32; c : out Complex_Number );

  -- DESCRIPTION :
  --   Using the given seed, generates a random complex number c,
  --   with real and imaginary parts both in [-1,+1].
  --   The seed is updated so the next call will give another c.

  function Random ( modulus : double_float ) return Complex_Number;

  -- DESCRIPTION :
  --   Generates a random complex number with a given modulus,
  --   so only the argument angle will be chosen at random.

  procedure Random_Complex_Number
              ( seed : in out integer32; modulus : in double_float;
                c : out Complex_Number );

  -- DESCRIPTION :
  --   Using the given seed, generates a random complex number c,
  --   with the given modulus, so only the argument is random.
  --   The seed is updated so the next call will give another c.

  function Random1 return Complex_Number;

  -- DESCRIPTION :
  --   Generates a random complex number with modulus one.

  procedure Random1_Complex_Number
              ( seed : in out integer32; c : out Complex_Number );

  -- DESCRIPTION :
  --   Given a seed, returns a random complex number c with modulus 1.
  --   The seed is updated so the next call returns a different c.

  function Random_Magnitude ( m : natural32 ) return double_float;
 
  -- DESCRIPTION :
  --   Returns a random double float in [10^(-m),10^(+m)].

  function Random_Magnitude ( m : natural32 ) return Complex_Number;

  -- DESCRIPTION :
  --   Returns a random complex number with modulus in [10^(-m),10^(+m)].

end Standard_Random_Numbers;
