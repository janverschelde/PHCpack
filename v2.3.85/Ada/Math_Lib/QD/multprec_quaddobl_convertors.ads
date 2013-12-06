with Multprec_Integer_Numbers;           use Multprec_Integer_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;

package Multprec_QuadDobl_Convertors is

-- DESCRIPTION :
--   Offers routine to convert a quad double in a multiprecision
--   floating number and to convert a multiprecision floating number
--   into a quad double.

  function to_floating_number ( d : quad_double ) return Floating_Number;

  -- DESCRIPTION :
  --   Returns a multprecision floating number with value equal
  --   to the quad double d.

  function to_quad_double ( i : Integer_Number ) return quad_double;

  -- DESCRIPTION :
  --   Returns the quad double that approximates i best.

  function to_quad_double ( f : Floating_Number ) return quad_double;

  -- DECRIPTION :
  --   Returns the quad double number that approximates
  --   the multiprecision number f best.

end Multprec_QuadDobl_Convertors;
