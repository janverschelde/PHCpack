with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Multprec_Natural_Numbers;
with Multprec_Natural64_Numbers;
with Multprec_Integer_Numbers;
with Multprec_Integer64_Numbers;
with Multprec_Floating_Numbers;
with Multprec_Floating64_Numbers;
with Multprec_Complex_Numbers;           use Multprec_Complex_Numbers;

package Multprec_Random_Numbers is

-- DESCRIPTION :
--   This package provides the generation of multi-precision natural
--   and integer number of a prescribed size.

  function Random ( size : natural32 )
                  return Multprec_Natural_Numbers.Natural_Number;
  function Random ( size : natural32 )
                  return Multprec_Natural64_Numbers.Natural_Number;

  -- DESCRIPTION :
  --   Returns a random natural number of the given size.

  function Random ( size : natural32 )
                  return Multprec_Integer_Numbers.Integer_Number;
  function Random ( size : natural32 )
                  return Multprec_Integer64_Numbers.Integer_Number;

  -- DESCRIPTION :
  --   Returns a random integer number of the given size.

  function Random ( size : natural32 )
                  return Multprec_Floating_Numbers.Floating_Number;
  function Random ( size : natural32 )
                  return Multprec_Floating64_Numbers.Floating_Number;

  -- DESCRIPTION :
  --   Returns a floating number of length size
  --   and the fraction in [0.1,1.0].

  function Random ( size : natural32; lower,upper : integer32 )
                  return Multprec_Floating_Numbers.Floating_Number;
  function Random ( size : natural32; lower,upper : integer64 )
                  return Multprec_Floating64_Numbers.Floating_Number;

  -- DESCRIPTION :
  --   Generates a random number of the given size, with fraction
  --   in [0.1,1.0] and the exponent between lower and upper.

  function Random ( size : natural32 ) return Complex_Number;

  -- DESCRIPTION :
  --   Returns a random complex number of the given size.

  function Random ( size : natural32; lower,upper : integer32 )
                  return Multprec_Complex_Numbers.Complex_Number;

  -- DESCRIPTION :
  --   Generates a random number of the given size,
  --   with exponents for real and imaginary parts between lower and upper. 
  --   Fractions lie in [0.1,1.0].

end Multprec_Random_Numbers;
