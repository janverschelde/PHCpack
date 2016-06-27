with DoblDobl_Dense_Series_Ring;
with DoblDobl_Series_Polynomials;
with Generic_Polynomial_Systems;

package DoblDobl_Series_Poly_Systems is
  new Generic_Polynomial_Systems(DoblDobl_Dense_Series_Ring,
                                 DoblDobl_Series_Polynomials);

-- DESCRIPTION :
--   Defines polynomial systems with polynomials with coefficients
--   as series of complex numbers, in double double precision.
