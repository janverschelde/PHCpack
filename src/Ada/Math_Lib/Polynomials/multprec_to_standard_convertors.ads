with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Multprec_Complex_Polynomials;
with Multprec_Complex_Poly_Systems;

package Multprec_to_Standard_Convertors is

-- DESCRIPTION :
--   Conversion from multi-precision complex into standard complex
--   for the coefficients of multivariate polynomials.
--   and setting of the size of the coefficients.

  function Convert ( p : Multprec_Complex_Polynomials.Poly )
                   return Standard_Complex_Polynomials.Poly;

  function Convert ( p : Multprec_Complex_Poly_Systems.Poly_Sys )
                   return Standard_Complex_Poly_Systems.Poly_Sys;

end Multprec_to_Standard_Convertors;
