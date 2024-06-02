with Standard_CSeries_Polynomials;
with Standard_CSeries_Poly_Systems;
with Double_Taylor_Homotopies;

package Taylor_Homotopy_Series is

-- DESCRIPTION :
--   Defines the conversion of a Taylor monomial homotopy into
--   a system with power series coefficients.

  function Make ( m : Double_Taylor_Homotopies.Taylor_Monomial )
                return Standard_CSeries_Polynomials.Term;
  function Make ( m : Double_Taylor_Homotopies.Link_to_Taylor_Monomial )
                return Standard_CSeries_Polynomials.Term;

  -- DESCRIPTION :
  --   Returns a monomial with series coefficient,
  --   given a Taylor monomial in m.

  function Make ( v : Double_Taylor_Homotopies.Taylor_Monomial_Vector )
                return Standard_CSeries_Polynomials.Poly;
  function Make ( v : Double_Taylor_Homotopies.Link_to_Taylor_Monomial_Vector )
                return Standard_CSeries_Polynomials.Poly;

  -- DESCRIPTION :
  --   Returns the polynomial with series coefficients from the monomials
  --   in the given vector v.

  function Make ( h : Double_Taylor_Homotopies.Taylor_Homotopy )
                return Standard_CSeries_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Returns the system of polynomials with series coefficients
  --   from the Taylor homotopy given in h.

end Taylor_Homotopy_Series;
