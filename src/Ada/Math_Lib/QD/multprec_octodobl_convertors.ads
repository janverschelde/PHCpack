with Multprec_Integer_Numbers;           use Multprec_Integer_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Octo_Double_Numbers;                use Octo_Double_Numbers;

package Multprec_OctoDobl_Convertors is

-- DESCRIPTION :
--   Offers routine to convert a octo double in a multiprecision
--   floating number and to convert a multiprecision floating number
--   into a octo double.

  function to_floating_number ( d : octo_double ) return Floating_Number;

  -- DESCRIPTION :
  --   Returns a multprecision floating number with value equal
  --   to the octo double d.

  function to_octo_double ( i : Integer_Number ) return octo_double;

  -- DESCRIPTION :
  --   Returns the octo double number that approximates i best.

  function to_octo_double ( f : Floating_Number ) return octo_double;

  -- DECRIPTION :
  --   Returns the octo double number that approximates
  --   the multiprecision number f best.

end Multprec_OctoDobl_Convertors;
