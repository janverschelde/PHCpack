with Multprec_Integer_Numbers;           use Multprec_Integer_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;

package Multprec_DoblDobl_Convertors is

-- DESCRIPTION :
--   Offers routine to convert a double double in a multiprecision
--   floating number and to convert a multiprecision floating number
--   into a double double.

  function to_floating_number ( d : double_double ) return Floating_Number;

  -- DESCRIPTION :
  --   Returns a multprecision floating number with value equal
  --   to the double double d.

  function to_double_double ( i : Integer_Number ) return double_double;

  -- DESCRIPTION :
  --   Returns the double double number that approximates i best.

  function to_double_double ( f : Floating_Number ) return double_double;

  -- DECRIPTION :
  --   Returns the double double number that approximates
  --   the multiprecision number f best.

end Multprec_DoblDobl_Convertors;
