with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Multprec_Integer_Numbers;           use Multprec_Integer_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;

package Standard_Complex_Exponentiation is

-- DESCRIPTION :
--   The operations in this package allow the calculation of powers of
--   standard complex numbers with very high exponents.

  procedure DivModTwoPi ( x : in double_float;
                          q : out integer32; r : out double_float );
  procedure DivModTwoPi ( x : in double_double;
                          q : out Integer_Number; r : out double_double );
  procedure DivModTwoPi ( x : in quad_double;
                          q : out Integer_Number; r : out quad_double );

  -- DESCRIPTION :
  --   Computes quotient q and remainder r after division of x by 2*Pi,
  --   so on return we have x = 2*Pi*q + r.
  --   If x > -2*Pi and x < 2*Pi, then q = 0 and r = x on return.

  procedure DivModTwoPi ( x : in Floating_Number; d : in natural32;
                          q : out Integer_Number; r : out Floating_Number );

  -- DESCRIPTION :
  --   Computes quotient q and remainder r after division of x by 2*Pi,
  --   using d decimal places of Pi.  On return we have x = 2*Pi + r.
  --   If x > -2*Pi and x < 2*Pi, then q = 0 and r = x on return.

  function Polar_Exponentiation_ModTwoPi
             ( x : Complex_Number; e : integer32 ) return Complex_Number;

  -- DESCRIPTION :
  --   Returns x^e via polar coordinates of x, computing the remainder
  --   of the exponent modulo 2*Pi.

  function Polar_Exponentiation_ModTwoPi_of_Unit
             ( x : Complex_Number; e : integer32 ) return Complex_Number;

  -- DESCRIPTION :
  --   Returns x^e via polar coordinates of x, where |x| = 1,
  --   computing the remainder of the exponent modulo 2*Pi.

  function Polar_Exponentiation_ModTwoPi_of_Unit
             ( x : Complex_Number; e : Integer_Number ) return Complex_Number;

  -- DESCRIPTION :
  --   Returns x^e via polar coordinates of x, where |x| = 1,
  --   computing the remainder of the exponent modulo 2*Pi.

end Standard_Complex_Exponentiation;
