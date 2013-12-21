with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with QuadDobl_Complex_Numbers;           use QuadDobl_Complex_Numbers;

package QuadDobl_Random_Numbers is

-- DESCRIPTION :
--   This package provides random number generators for 
--   quad doubles and complex quad doubles.
--   The seed is set in the package Standard_Random_Numbers.

  function Random return quad_double;

  -- DESCRIPTION :
  --   Returns a random quad double in [-1,+1].

  function Random_Magnitude ( m : natural32 ) return quad_double;

  -- DESCRIPTION :
  --   Returns a random quad double with absolute value
  --   in [10^(-m),10^(+m)].

  function Random return Complex_Number;

  -- DESCRIPTION :
  --   Returns a random complex quad double with
  --   real and imaginary parts in [-1,+1].

  function Random1 return Complex_Number;

  -- DESCRIPTION :
  --   Returns a random complex quad double
  --   with modulus equal to one.

  function Random_Magnitude ( m : natural32 ) return Complex_Number;

  -- DESCRIPTION :
  --   Returns a random complex quad double of modulus in
  --   the interval [10^(-m),10^(+m)].

end QuadDobl_Random_Numbers;
