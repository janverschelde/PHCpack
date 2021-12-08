with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Hexa_Double_Numbers;                use Hexa_Double_Numbers;
with HexaDobl_Complex_Numbers;           use HexaDobl_Complex_Numbers;

package HexaDobl_Random_Numbers is

-- DESCRIPTION :
--   This package provides random number generators for 
--   hexa doubles and complex hexa doubles.
--   The seed is set in the package Standard_Random_Numbers.
--   Alternatively, the user can manage the seed for independent
--   and/or reproducible sequences of random numbers.

  function Random return hexa_double;

  -- DESCRIPTION :
  --   Returns a random hexa double in [-1,+1].

  procedure Random_Hexa_Double
              ( seed : in out integer32; f : out hexa_double );

  -- DESCRIPTION :
  --   Given a seed, returns a random hexa double f in [-1,+1].
  --   The seed is updated so the next call returns a different f.
 
  function Random return Complex_Number;

  -- DESCRIPTION :
  --   Returns a random complex hexa double,
  --   with real and imaginary parts in [-1,+1].

  procedure Random_Complex_Number
              ( seed : in out integer32; c : out Complex_Number );

  -- DESCRIPTION :
  --   Given a seed, returns a random complex number c,
  --   with real and imaginary parts both in [-1,+1].
  --   The seed is updated so the next call returns a different c.

  function Random1 return Complex_Number;

  -- DESCRIPTION :
  --   Returns a random complex hexa double with modulus equal to one.

  procedure Random1_Complex_Number
              ( seed : in out integer32; c : out Complex_Number );

  -- DESCRIPTION :
  --   Given a seed, returns a random complex number c,
  --   with modulus equal to one.
  --   The seed is updated so the next call returns a different c.

end HexaDobl_Random_Numbers; 
