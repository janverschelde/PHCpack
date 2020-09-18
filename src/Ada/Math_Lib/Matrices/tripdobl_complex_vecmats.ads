with Generic_VecMats;
with TripDobl_Complex_Ring;
with TripDobl_Complex_Vectors;
with TripDobl_Complex_Matrices;

package TripDobl_Complex_VecMats is 
  new Generic_VecMats(TripDobl_Complex_Ring,
                      TripDobl_Complex_Vectors,
                      TripDobl_Complex_Matrices);

-- DESCRIPTION :
--   Defines vectors of matrices over the ring of complex numbers
--   with real and imaginary triple double floating-point numbers.
