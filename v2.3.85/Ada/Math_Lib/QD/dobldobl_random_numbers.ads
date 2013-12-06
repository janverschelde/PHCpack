with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with DoblDobl_Complex_Numbers;           use DoblDobl_Complex_Numbers;

package DoblDobl_Random_Numbers is

-- DESCRIPTION :
--   This package provides random number generators for 
--   double doubles and complex double doubles.
--   The seed is set in the package Standard_Random_Numbers.

  function Random return double_double;

  -- DESCRIPTION :
  --   Returns a random double double in [-1,+1].
 
  function Random_Magnitude ( m : natural32 ) return double_double;

  -- DESCRIPTION :
  --   Returns a random double double with absolute value
  --   in [10^(-m),10^(+m)].

  function Random return Complex_Number;

  -- DESCRIPTION :
  --   Returns a random complex double double with
  --   real and imaginary parts in [-1,+1].

  function Random1 return Complex_Number;

  -- DESCRIPTION :
  --   Returns a random complex double double
  --   with modulus equal to one.
 
  function Random_Magnitude ( m : natural32 ) return Complex_Number;

  -- DESCRIPTION :
  --   Returns a random complex double double of modulus in
  --   the interval [10^(-m),10^(+m)].

end DoblDobl_Random_Numbers;
