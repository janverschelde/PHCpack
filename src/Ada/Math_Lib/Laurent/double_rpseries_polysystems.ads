with Generic_Laur_Poly_Systems;
with Double_rpSeries_Ring;
with Double_rpSeries_Polynomials;

package Double_rpSeries_PolySystems is 
  new Generic_Laur_Poly_Systems(Double_rpSeries_Ring,
                                Double_rpSeries_Polynomials);

-- DESCRIPTION :
--   Defines the Laurent polynomial systems over the ring of series 
--   with real exponents in double precision.
