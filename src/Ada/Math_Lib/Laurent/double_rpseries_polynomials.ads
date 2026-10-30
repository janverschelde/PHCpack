with Generic_Laurent_Polynomials;
with Double_rpSeries_Ring;

package Double_rpSeries_Polynomials is 
  new Generic_Laurent_Polynomials(Double_rpSeries_Ring);

-- DESCRIPTION :
--   Defines the Laurent polynomials over the ring of series 
--   with real exponents in double precision.
