with Double_rpSeries_Ring;
with Double_rpSeries_Ring.FField;
with Double_rpSeries_Vectors;
with Double_rpSeries_Polynomials;
with Generic_Laur_Poly_Functions;

package Double_rpSeries_PolyFunctions is
  new Generic_Laur_Poly_Functions(Double_rpSeries_Ring,
                                  Double_rpSeries_Ring.FField,
                                  Double_rpSeries_Vectors,
                                  Double_rpSeries_Polynomials);

-- DESCRIPTION :
--   Defines Laurent polynomial functions for series with real exponents
--   and complex coefficients in double precision.
