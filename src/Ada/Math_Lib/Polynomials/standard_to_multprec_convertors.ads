with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Polynomials;
with Standard_Complex_Laurentials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with Multprec_Floating_Polynomials;
with Multprec_Floating_Poly_Systems;
with Multprec_Complex_Laurentials;
with Multprec_Complex_Polynomials;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Laur_Systems;

package Standard_to_Multprec_Convertors is

-- DESCRIPTION :
--   Conversion from standard complex into multi-precision complex
--   and setting of the size of the coefficients.

  function Convert ( p : Standard_Complex_Polynomials.Poly )
                   return Multprec_Complex_Polynomials.Poly;
  function Convert ( p : Standard_Complex_Laurentials.Poly )
                   return Multprec_Complex_Laurentials.Poly;

  -- DESCRIPTION :
  --   Converts the polynomial p with standard complex floating-point
  --   coefficients into a corresponding multiprecision polynomial.

  function Convert ( p : Standard_Complex_Poly_Systems.Poly_Sys )
                   return Multprec_Complex_Poly_Systems.Poly_Sys;
  function Convert ( p : Standard_Complex_Laur_Systems.Laur_Sys )
                   return Multprec_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION :
  --   Returns a system of converted polynomials.

  procedure Set_Size ( p : in out Multprec_Floating_Polynomials.Poly;
                       size : in natural32 );
  procedure Set_Size ( p : in out Multprec_Complex_Polynomials.Poly;
                       size : in natural32 );
  procedure Set_Size ( p : in out Multprec_Complex_Laurentials.Poly;
                       size : in natural32 );
  procedure Set_Size ( p : in out Multprec_Floating_Poly_Systems.Poly_Sys;
                       size : in natural32 );
  procedure Set_Size ( p : in out Multprec_Complex_Poly_Systems.Poly_Sys;
                       size : in natural32 );
  procedure Set_Size ( p : in out Multprec_Complex_Laur_Systems.Laur_Sys;
                       size : in natural32 );

  -- DESCRIPTION :
  --   Sets the size of all coefficients in p to the given size.

end Standard_to_Multprec_Convertors;
