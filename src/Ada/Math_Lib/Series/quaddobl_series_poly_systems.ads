with QuadDobl_Dense_Series_Ring;
with QuadDobl_Series_Polynomials;
with Generic_Polynomial_Systems;

package QuadDobl_Series_Poly_Systems is
  new Generic_Polynomial_Systems(QuadDobl_Dense_Series_Ring,
                                 QuadDobl_Series_Polynomials);

-- DESCRIPTION :
--   Defines polynomial systems with polynomials with coefficients
--   as series of complex numbers, in quad double precision.
