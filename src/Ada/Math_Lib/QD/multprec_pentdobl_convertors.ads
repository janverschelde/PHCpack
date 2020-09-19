with Multprec_Integer_Numbers;           use Multprec_Integer_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Penta_Double_Numbers;               use Penta_Double_Numbers;

package Multprec_PentDobl_Convertors is

-- DESCRIPTION :
--   Offers routine to convert a penta double in a multiprecision
--   floating number and to convert a multiprecision floating number
--   into a penta double.

  function to_floating_number ( d : penta_double ) return Floating_Number;

  -- DESCRIPTION :
  --   Returns a multprecision floating number with value equal
  --   to the penta double d.

  function to_penta_double ( i : Integer_Number ) return penta_double;

  -- DESCRIPTION :
  --   Returns the penta double number that approximates i best.

  function to_penta_double ( f : Floating_Number ) return penta_double;

  -- DECRIPTION :
  --   Returns the penta double number that approximates
  --   the multiprecision number f best.

end Multprec_PentDobl_Convertors;
