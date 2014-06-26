with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with DoblDobl_Complex_Numbers;           use DoblDobl_Complex_Numbers;
with Multprec_Integer_Numbers;           use Multprec_Integer_Numbers;

package DoblDobl_Complex_Exponentiation is

-- DESCRIPTION :
--   The operations in this package allow the calculation of powers of
--   standard complex numbers with very high exponents.

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

end DoblDobl_Complex_Exponentiation;
