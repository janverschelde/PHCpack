with Multprec_Integer_Numbers;           use Multprec_Integer_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Triple_Double_Numbers;              use Triple_Double_Numbers;

package Multprec_TripDobl_Convertors is

-- DESCRIPTION :
--   Offers routine to convert a triple double in a multiprecision
--   floating number and to convert a multiprecision floating number
--   into a triple double.

  function to_floating_number ( d : triple_double ) return Floating_Number;

  -- DESCRIPTION :
  --   Returns a multprecision floating number with value equal
  --   to the triple double d.

  function to_triple_double ( i : Integer_Number ) return triple_double;

  -- DESCRIPTION :
  --   Returns the triple double number that approximates i best.

  function to_triple_double ( f : Floating_Number ) return triple_double;

  -- DECRIPTION :
  --   Returns the triple double number that approximates
  --   the multiprecision number f best.

end Multprec_TripDobl_Convertors;
