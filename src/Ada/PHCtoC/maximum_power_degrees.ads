with Standard_Integer_Numbers;         use Standard_Integer_Numbers;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;

package Maximum_Power_Degrees is

-- DESCRIPTION :
--   Returns the maximum power of any variables in the degrees
--   of the monomials in a polynomial system.

  function Maximum_Power
             ( t : Standard_Complex_Polynomials.Term ) return integer32;
  function Maximum_Power
             ( t : DoblDobl_Complex_Polynomials.Term ) return integer32;
  function Maximum_Power
             ( t : QuadDobl_Complex_Polynomials.Term ) return integer32;

  -- DESCRIPTION :
  --   Returns the largest power in the degrees of t,
  --   over all variables in t.
  --   Returns -1 if the the term is empty.

  function Maximum_Power
             ( p : Standard_Complex_Polynomials.Poly ) return integer32;
  function Maximum_Power
             ( p : DoblDobl_Complex_Polynomials.Poly ) return integer32;
  function Maximum_Power
             ( p : QuadDobl_Complex_Polynomials.Poly ) return integer32;

  -- DESCRIPTION :
  --   Returns the largest power in the degrees of p.

  function Maximum_Power
             ( p : Standard_Complex_Poly_Systems.Poly_Sys ) return integer32;
  function Maximum_Power
             ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys ) return integer32;
  function Maximum_Power
             ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys ) return integer32;

  -- DESCRIPTION :
  --   Returns the largest power in the degrees of p.

end Maximum_Power_Degrees;
