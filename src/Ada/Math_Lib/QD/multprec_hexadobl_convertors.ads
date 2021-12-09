with Multprec_Integer_Numbers;           use Multprec_Integer_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Hexa_Double_Numbers;                use Hexa_Double_Numbers;

package Multprec_HexaDobl_Convertors is

-- DESCRIPTION :
--   Offers routine to convert a hexa double in a multiprecision
--   floating number and to convert a multiprecision floating number
--   into a hexa double.

  function to_floating_number ( d : hexa_double ) return Floating_Number;

  -- DESCRIPTION :
  --   Returns a multprecision floating number with value equal
  --   to the hexa double d.

  function to_hexa_double ( i : Integer_Number ) return hexa_double;

  -- DESCRIPTION :
  --   Returns the hexa double number that approximates i best.

  function to_hexa_double ( f : Floating_Number ) return hexa_double;

  -- DECRIPTION :
  --   Returns the hexa double number that approximates
  --   the multiprecision number f best.

end Multprec_HexaDobl_Convertors;
