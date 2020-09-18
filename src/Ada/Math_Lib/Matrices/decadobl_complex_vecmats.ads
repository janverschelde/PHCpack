with Generic_VecMats;
with DecaDobl_Complex_Ring;
with DecaDobl_Complex_Vectors;
with DecaDobl_Complex_Matrices;

package DecaDobl_Complex_VecMats is 
  new Generic_VecMats(DecaDobl_Complex_Ring,
                      DecaDobl_Complex_Vectors,
                      DecaDobl_Complex_Matrices);

-- DESCRIPTION :
--   Defines vectors of matrices over the ring of complex numbers
--   with real and imaginary deca double floating-point numbers.
