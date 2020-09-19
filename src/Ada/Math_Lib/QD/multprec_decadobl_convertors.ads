with Multprec_Integer_Numbers;           use Multprec_Integer_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Deca_Double_Numbers;                use Deca_Double_Numbers;

package Multprec_DecaDobl_Convertors is

-- DESCRIPTION :
--   Offers routine to convert a deca double in a multiprecision
--   floating number and to convert a multiprecision floating number
--   into a deca double.

  function to_floating_number ( d : deca_double ) return Floating_Number;

  -- DESCRIPTION :
  --   Returns a multprecision floating number with value equal
  --   to the deca double d.

  function to_deca_double ( i : Integer_Number ) return deca_double;

  -- DESCRIPTION :
  --   Returns the deca double number that approximates i best.

  function to_deca_double ( f : Floating_Number ) return deca_double;

  -- DECRIPTION :
  --   Returns the deca double number that approximates
  --   the multiprecision number f best.

end Multprec_DecaDobl_Convertors;
