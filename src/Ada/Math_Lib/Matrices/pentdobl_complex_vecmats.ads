with Generic_VecMats;
with PentDobl_Complex_Ring;
with PentDobl_Complex_Vectors;
with PentDobl_Complex_Matrices;

package PentDobl_Complex_VecMats is 
  new Generic_VecMats(PentDobl_Complex_Ring,
                      PentDobl_Complex_Vectors,
                      PentDobl_Complex_Matrices);

-- DESCRIPTION :
--   Defines vectors of matrices over the ring of complex numbers
--   with real and imaginary penta double floating-point numbers.
